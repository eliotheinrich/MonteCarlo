#include "Spin2DModel.h"
#include <math.h>
#include <fstream>

#define DEFAULT_CLUSTER_UPDATE true

#define DEFAULT_SAMPLE_ENERGY true
#define DEFAULT_SAMPLE_MAGNETIZATION true
#define DEFAULT_SAMPLE_HELICITY false

#define DEFAULT_BOUNDARY_CONDITION "periodic"

#define GHOST -1

using namespace dataframe::utils;

Spin2DModel::Spin2DModel(dataframe::Params &params, uint32_t num_threads) : MonteCarloSimulator(params, num_threads) {
  cluster_update = dataframe::utils::get<int>(params, "cluster_update", DEFAULT_CLUSTER_UPDATE);

  sample_energy = dataframe::utils::get<int>(params, "sample_energy", DEFAULT_SAMPLE_ENERGY);
  sample_magnetization = dataframe::utils::get<int>(params, "sample_magnetization", DEFAULT_SAMPLE_MAGNETIZATION);
  sample_helicity = dataframe::utils::get<int>(params, "sample_helicity", DEFAULT_SAMPLE_HELICITY);


  bcx = parse_boundary_condition(dataframe::utils::get<std::string>(params, "bcx", DEFAULT_BOUNDARY_CONDITION));
  bcy = parse_boundary_condition(dataframe::utils::get<std::string>(params, "bcy", DEFAULT_BOUNDARY_CONDITION));
  bcz = parse_boundary_condition(dataframe::utils::get<std::string>(params, "bcz", DEFAULT_BOUNDARY_CONDITION));
}

void Spin2DModel::init(uint32_t sl, uint32_t N1, uint32_t N2=-1, uint32_t N3=-1) {
  this->sl = sl;
  this->N1 = N1;

  if (N2 == -1) { 
    this->N2 = N1; 
  } else { 
    this->N2 = N2; 
  }

  if (N3 == -1) { 
    this->N3 = N1; 
  } else { 
    this->N3 = N3; 
  }

  V = N1*N2*N3*sl;

  spins = std::vector<Eigen::Vector2d>(V);
  if (cluster_update) {
    neighbors = std::vector<std::vector<Bond>>(V+1, std::vector<Bond>(0));
  } else {
    neighbors = std::vector<std::vector<Bond>>(V, std::vector<Bond>(0));
  }

  for (uint32_t n = 0; n < bonds.size(); n++) {
    auto b = bonds[n];
    auto bond_filter = bond_filters[n];

    for (uint32_t i = 0; i < V; i++) {
      Eigen::Vector4i idx = tensor_idx(i);

      uint nx = idx[0] + b.d1;
      uint ny = idx[1] + b.d2;
      uint nz = idx[2] + b.d3;
      uint ns = idx[3] + b.ds;

      if (bcx == BoundaryCondition::Open) {
        if (nx < 0 || nx >= N1) {
          continue;
        }
      } else if (bcx == BoundaryCondition::Periodic) {
        nx = mod(nx, N1);
      }

      if (bcy == BoundaryCondition::Open) {
        if (ny < 0 || ny >= N2) {
          continue;
        }
      } else if (bcx == BoundaryCondition::Periodic) {
        ny = mod(ny, N2);
      }

      if (bcz == BoundaryCondition::Open) {
        if (nz < 0 || nz >= N3) {
          continue;
        }
      } else if (bcz == BoundaryCondition::Periodic) {
        nz = mod(nz, N3);
      }

      uint j = flat_idx(nx, ny, nz, ns);
      if (!bond_filter(i, j)) {
        continue;
      }

      neighbors[i].push_back(Bond{j, n});
    }
  }

  if (cluster_update) {
    // Connect every site to the ghost 
    for (uint32_t i = 0; i < V; i++) {
      neighbors[V].push_back(Bond{i, GHOST});
      neighbors[i].push_back(Bond{V, GHOST});
    }

    s0 = Eigen::Matrix2d::Identity();
  }  

  randomize_spins();

  mut.i = 0;
  sigma = 0.25;
}

void Spin2DModel::randomize_spins() {
  double p;
  Eigen::Vector2d v;

  for (uint32_t i = 0; i < V; i++) {
    p = 2*PI*randf();
    spins[i] << std::cos(p), std::sin(p);
  }
}

void Spin2DModel::add_bond(
  int d1, 
  int d2, 
  int d3, 
  int ds, 
  Eigen::Vector3d v, 
  std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)> bondfunc,
  std::function<bool(uint32_t, uint32_t)> bond_filter
) {
  Spin2DBond b{d1, d2, d3, ds, v, bondfunc};
  bonds.push_back(b);
  bond_filters.push_back(bond_filter);

  Eigen::Matrix2d R;
  R << std::cos(v[0]*alpha), -std::sin(v[0]*alpha),
       std::sin(v[0]*alpha),  std::cos(v[0]*alpha);
  R1s.push_back(R);
  R2s.push_back(R*R);
  R3s.push_back(R*R*R);
}

std::vector<double> Spin2DModel::twist_stiffness() const {
  // Returns the first and second derivative in response to a phase twist
  double E0 = 0.;
  double E1 = 0.;
  double E2 = 0.;
  double E3 = 0.;
  double Em1 = 0.;
  double Em2 = 0.;
  double Em3 = 0.;


  Eigen::Vector2d S1;
  Eigen::Vector2d S2;
  for (uint32_t i = 0; i < V; i++) {
    for (auto const &[j, b] : neighbors[i]) {
      S1 = spins[i];
      S2 = spins[j];

      E0 += bonds[b].bondfunc(S1, S2);

      E1 += bonds[b].bondfunc(S1, R1s[b]*S2);
      Em1 += bonds[b].bondfunc(S1, R1s[b].transpose()*S2);

      E2 += bonds[b].bondfunc(S1, R2s[b]*S2);
      Em2 += bonds[b].bondfunc(S1, R2s[b].transpose()*S2);

      E3 += bonds[b].bondfunc(S1, R3s[b]*S2);
      Em3 += bonds[b].bondfunc(S1, R3s[b].transpose()*S2);
    }
  }

  // Compute derivates from finite difference
  double d1E = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha/2.;
  double d2E = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/pow(alpha, 2)/2.;
  double d3E = (1./8.*Em3 - 1.*Em2 + 13./8.*Em1 - 13./8.*E1 + 1.*E2 - 1./8.*E3)/pow(alpha, 3)/2.;
  double d4E = (-1./6.*Em3 + 2.*Em2 - 13./2.*Em1 + 28./3.*E0 - 13./2.*E1 + 2.*E2 - 1./6.*E3)/pow(alpha, 4)/2.;

  return std::vector<double>{d4E, d3E, d1E, pow(d2E,2), d2E, d1E*d3E, pow(d1E,2)*d2E, 
    d1E*d2E, pow(d1E,2), pow(d1E,4), pow(d1E,3)};
}

Eigen::Vector2d Spin2DModel::get_magnetization() const {
  Eigen::Vector2d M = Eigen::Vector2d::Constant(0);
  for (uint32_t i = 0; i < V; i++) {
    M += spins[i];
  }

  if (cluster_update) {
    return s0.transpose()*M/V;
  } else {
    return M/V;
  }
}

void Spin2DModel::generate_mutation() {
  if (cluster_update) {
    cluster_mutation();
  } else {
    mut.i = rand() % V;
    metropolis_mutation();
  }
}

void Spin2DModel::cluster_mutation() {
  s.clear();

  std::stack<uint32_t> c;
  uint32_t m = rand() % (V + 1);
  Eigen::Matrix2d s0_new;
  c.push(m);

  double p = (double) 2*PI*rand()/RAND_MAX;
  Eigen::Vector2d ax; ax << std::cos(p), std::sin(p);
  Eigen::Matrix2d R = Eigen::Matrix2d::Identity() - 2*ax*ax.transpose()/std::pow(ax.norm(), 2);

  Eigen::Vector2d s_new;
  while (!c.empty()) {
    m = c.top();
    c.pop();

    if (!s.count(m)) {
      s.insert(m);
      bool is_ghost = (m == V);
      if (is_ghost) { // Site is ghost
        s0_new = R*s0;
      } else {
        s_new = R*spins[m];
      }

      for (auto const &[j, b] : neighbors[m]) {
        if (!s.count(j)) {
          bool neighbor_is_ghost = (j == V);
          
          double dE;
          if (neighbor_is_ghost) {
            dE = onsite_func(s0.inverse()*s_new) - onsite_func(s0.inverse()*spins[m]);
          } else if (is_ghost) {
            dE = onsite_func(s0_new.inverse()*spins[j]) - onsite_func(s0.inverse()*spins[j]);
          } else { // Normal bond
            dE = bonds[b].bondfunc(spins[j], s_new) - bonds[b].bondfunc(spins[j], spins[m]);
          }

          if (randf() < 1. - std::exp(-dE/temperature)) {
            c.push(j);
          }
        }
      }

      if (is_ghost) {
        s0 = s0_new;
      } else {
        spins[m] = s_new;
      }
    }
  }
}

void Spin2DModel::metropolis_mutation() {
  double dp = sigma*(randf() - 0.5)*2.*PI;
  Eigen::Vector2d S1 = spins[mut.i];
  Eigen::Vector2d S2; S2 << std::cos(dp)*S1[0] - std::sin(dp)*S1[1],
                            std::cos(dp)*S1[1] + std::sin(dp)*S1[0];

  // Store mutation for consideration
  mut.dS = S2 - S1;
}


void Spin2DModel::accept_mutation() {
  return;
}

void Spin2DModel::reject_mutation() {
  if (!cluster_update) {
    spins[mut.i] = spins[mut.i] - mut.dS;
  }
}

double Spin2DModel::onsite_energy(int i) const {
  if (cluster_update) {
    return onsite_func(s0.transpose()*spins[i]);
  } else {
    return onsite_func(spins[i]);
  }
}

double Spin2DModel::bond_energy(int i) const {
  double E = 0.;
  for (auto const &[j, b] : neighbors[i]) {
    if (b == GHOST) {
      continue;
    }

    E += 0.5*bonds[b].bondfunc(spins[i], spins[j]);
  }

  return E;
}

double Spin2DModel::energy() const {
  double E = 0;

  for (uint32_t i = 0; i < V; i++) {
    E += onsite_energy(i);
    E += bond_energy(i);
  }

  return E;
}

double Spin2DModel::energy_change() {
  if (cluster_update) {
    return -1.;
  }

  double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
  spins[mut.i] = spins[mut.i] + mut.dS;
  double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

  return E2 - E1;
}

void Spin2DModel::save_spins(const std::string& filename) {
  std::ofstream output_file;
  output_file.open(filename);
  output_file << N1 << "\t" << N2 << "\t" << N3 << "\t" << sl << "\n";

  Eigen::Vector2d S;
  for (uint32_t i = 0; i < V; i++) {
    if (cluster_update) {
      S = s0.transpose()*spins[i];
    } else {
      S = spins[i];
    }

    output_file << S[0] << "\t" << S[1];
    if (i < V-1) { 
      output_file << "\t"; 
    }
  }

  output_file.close();
}

void Spin2DModel::add_magnetization_samples(dataframe::data_t &samples) const {
  auto m = get_magnetization();
  samples.emplace("mx", m[0]);
  samples.emplace("my", m[1]);
  samples.emplace("magnetization", m.norm());
}

void Spin2DModel::add_helicity_samples(dataframe::data_t &samples) const {
  std::vector<double> twistd = twist_stiffness();
  samples.emplace("d1E", twistd[2]);
  samples.emplace("d2E", twistd[4]);
}

dataframe::data_t Spin2DModel::take_samples() {
  dataframe::data_t samples;

  if (sample_energy) {
    samples.emplace("energy", energy());
  }

  if (sample_magnetization) {
    add_magnetization_samples(samples);
  }

  if (sample_helicity) {
    add_helicity_samples(samples);
  }

  return samples;
}

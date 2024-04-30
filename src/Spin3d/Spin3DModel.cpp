#include "Spin3DModel.h"
#include <fstream>
#include <stack>

#define DEFAULT_CLUSTER_UPDATE true

#define DEFAULT_SAMPLE_ENERGY true
#define DEFAULT_SAMPLE_MAGNETIZATION true
#define DEFAULT_SAMPLE_HELICITY false

#define DEFAULT_BOUNDARY_CONDITION "periodic"


#define GHOST -1

GaussianDist::GaussianDist(double mean, double std) {
  rd.seed(rand());
  gen = std::default_random_engine(rd());
  dist = std::normal_distribution<>(mean, std);
}

double GaussianDist::sample() {
  return dist(gen);
}

Spin3DModel::Spin3DModel(dataframe::Params &params, uint32_t num_threads) : MonteCarloSimulator(params, num_threads), nsteps(0), accepted(0) {
  cluster_update = dataframe::utils::get<int>(params, "cluster_update", DEFAULT_CLUSTER_UPDATE);

  sample_helicity = dataframe::utils::get<int>(params, "sample_helicity", DEFAULT_SAMPLE_HELICITY);
  sample_magnetization = dataframe::utils::get<int>(params, "sample_magnetization", DEFAULT_SAMPLE_MAGNETIZATION);
  sample_energy = dataframe::utils::get<int>(params, "sample_energy", DEFAULT_SAMPLE_ENERGY);

  bcx = parse_boundary_condition(dataframe::utils::get<std::string>(params, "bcx", DEFAULT_BOUNDARY_CONDITION));
  bcy = parse_boundary_condition(dataframe::utils::get<std::string>(params, "bcy", DEFAULT_BOUNDARY_CONDITION));
  bcz = parse_boundary_condition(dataframe::utils::get<std::string>(params, "bcz", DEFAULT_BOUNDARY_CONDITION));
}

void Spin3DModel::init(uint32_t sl, uint32_t N1, uint32_t N2=-1, uint32_t N3=-1) {
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

  spins = std::vector<Eigen::Vector3d>(V);
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

    s0 = Eigen::Matrix3d::Identity();
  }

  randomize_spins();

  acceptance = 0.5;
  sigma = 0.25;

  dist = GaussianDist(0., 1.0);

  mut.i = 0;
}

void Spin3DModel::randomize_spins() {
  for (uint32_t i = 0; i < V; i++) {
    spins[i] = Eigen::Vector3d::Random(3).normalized();
  }
}

void Spin3DModel::add_bond(
  int d1, 
  int d2, 
  int d3, 
  int ds, 
  Eigen::Vector3d v, 
  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc, 
  std::function<bool(uint32_t, uint32_t)> bond_filter
) {
  Spin3DBond b{d1, d2, d3, ds, v, bondfunc};
  bonds.push_back(b);
  bond_filters.push_back(bond_filter);

  double f = v[0]*alpha;
  Eigen::Matrix3d R;
  R << std::cos(f), -std::sin(f), 0,
       std::sin(f),  std::cos(f), 0.,
       0.,           0.,          1.;
  R1s.push_back(R);
  R2s.push_back(R*R);
  R3s.push_back(R*R*R);
}

std::vector<double> Spin3DModel::twist_terms(std::vector<double> dE) {
  return std::vector<double>{dE[3], dE[2], dE[0], pow(dE[1],2), dE[1], dE[0]*dE[2], std::pow(dE[0],2)*dE[1], 
    dE[0]*dE[1], std::pow(dE[0],2), std::pow(dE[0],4), std::pow(dE[0],3)};
}

std::vector<double> Spin3DModel::twist_derivatives(uint32_t i) const {
  double E0 = 0.;
  double E1 = 0.;
  double E2 = 0.;
  double E3 = 0.;
  double Em1 = 0.;
  double Em2 = 0.;
  double Em3 = 0.;

  Eigen::Vector3d S1 = get_spin(i);
  Eigen::Vector3d S2;
  for (auto const &[j, b] : neighbors[i]) {
    if (b == GHOST) {
      continue;
    }

    S2 = get_spin(j);

    E0 += bonds[b].bondfunc(S1, S2);

    E1 += bonds[b].bondfunc(S1, R1s[b]*S2);
    Em1 += bonds[b].bondfunc(S1, R1s[b].transpose()*S2);

    E2 += bonds[b].bondfunc(S1, R2s[b]*S2);
    Em2 += bonds[b].bondfunc(S1, R2s[b].transpose()*S2);

    E3 += bonds[b].bondfunc(S1, R3s[b]*S2);
    Em3 += bonds[b].bondfunc(S1, R3s[b].transpose()*S2);
  }

  double d1E = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha/2.;
  double d2E = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/pow(alpha, 2)/2.;
  double d3E = (1./8.*Em3 - 1.*Em2 + 13./8.*Em1 - 13./8.*E1 + 1.*E2 - 1./8.*E3)/pow(alpha, 3)/2.;
  double d4E = (-1./6.*Em3 + 2.*Em2 - 13./2.*Em1 + 28./3.*E0 - 13./2.*E1 + 2.*E2 - 1./6.*E3)/pow(alpha, 4)/2.;

  return std::vector<double>{d1E, d2E, d3E, d4E};
}

std::vector<double> Spin3DModel::twist_derivatives() const {
  std::vector<double> twist = std::vector<double>(4, 0);
  std::vector<double> twist_i = std::vector<double>(4, 0);

  for (uint32_t i = 0; i < V; i++) {
    twist_i = twist_derivatives(i);
    for (uint32_t j = 0; j < 4; j++) {
      twist[j] += twist_i[j];
    }
  }

  return twist;
}

Eigen::Vector3d Spin3DModel::get_magnetization() const {
  Eigen::Vector3d M = Eigen::Vector3d::Constant(0);
  for (uint32_t i = 0; i < V; i++) {
    M += get_spin(i);
  }

  return M/V;
}

std::vector<double> Spin3DModel::correlation_function(uint32_t i, uint32_t a = 2, uint32_t b = 2) const {
  std::vector<double> Cij = std::vector<double>(V); 

  Eigen::Vector4i idxs = tensor_idx(i);
  int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];
  for (uint32_t n1 = 0; n1 < N1; n1++) {
    for (uint32_t n2 = 0; n2 < N2; n2++) {
      for (uint32_t n3 = 0; n3 < N3; n3++) {
        for (uint32_t s = 0; s < sl; s++) {
          uint32_t j = flat_idx(n1, n2, n3, s);
          Cij[j] = spins[flat_idx((m1 + n1)%N1, 
                                  (m2 + n2)%N2, 
                                  (m3 + n3)%N3, 
                                  (s + k)%sl)][a]*spins[i][b];
        }
      }
    }
  }
  return Cij;
}

std::vector<double> Spin3DModel::full_correlation_function(uint32_t i) const {
  std::vector<double> Cij = std::vector<double>(V); 

  Eigen::Vector4i idxs = tensor_idx(i);
  int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];
  for (uint32_t n1 = 0; n1 < N1; n1++) {
    for (uint32_t n2 = 0; n2 < N2; n2++) {
      for (uint32_t n3 = 0; n3 < N3; n3++) {
        for (uint32_t s = 0; s < sl; s++) {
          uint32_t j = flat_idx(n1, n2, n3, s);
          Cij[j] = spins[flat_idx((m1 + n1)%N1, 
                                  (m2 + n2)%N2, 
                                  (m3 + n3)%N3, 
                                  (s + k)%sl)].dot(spins[i]);
        }
      }
    }
  }
  return Cij;
}

double Spin3DModel::skyrmion_density(uint32_t i) const {
  Eigen::Vector3d dSdX; dSdX << 0., 0., 0.;
  Eigen::Vector3d dSdY; dSdY << 0., 0., 0.;
  for (auto const &[j, b] : neighbors[i]) {
    if (b == GHOST) {
      continue;
    }

    if (bonds[b].v[0] != 0.) {
      dSdX += bonds[b].v[0]*(spins[j] - spins[i]);
    }

    if (bonds[b].v[1] != 0.) {
      dSdX += bonds[b].v[1]*(spins[j] - spins[i]);
    }
  }
  dSdX = dSdX/bonds.size();
  dSdY = dSdY/bonds.size();

  return spins[i].dot(dSdX.cross(dSdY));
}

std::vector<double> Spin3DModel::skyrmion_correlation_function(uint32_t i) const {
  std::vector<double> Cij = std::vector<double>(V); 

  int j;

  Eigen::Vector4i idxs = tensor_idx(i);
  int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];

  double Si = skyrmion_density(i);
  for (uint32_t n1 = 0; n1 < N1; n1++) {
    for (uint32_t n2 = 0; n2 < N2; n2++) {
      for (uint32_t n3 = 0; n3 < N3; n3++) {
        for (uint32_t s = 0; s < sl; s++) {
          j = flat_idx(n1, n2, n3, s);
          Cij[j] = skyrmion_density(flat_idx((m1 + n1)%N1, 
                                             (m2 + n2)%N2, 
                                             (m3 + n3)%N3, 
                                             (s + k)%sl))*Si;
        }
      }
    }
  }
  return Cij;
}

void Spin3DModel::cluster_mutation() {
  s.clear();

  std::stack<uint32_t> c;
  uint32_t m = rand() % V;
  c.push(m);

  Eigen::Vector3d ax; ax << dist.sample(), dist.sample(), dist.sample();
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity() - 2*ax*ax.transpose()/std::pow(ax.norm(), 2);

  Eigen::Matrix3d s0_new;
  Eigen::Vector3d s_new;
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
        bool neighbor_is_ghost = (b == GHOST);
        if (!s.count(j)) {
          double dE;
          if (neighbor_is_ghost) {
            dE = onsite_func(s0.transpose()*s_new) - onsite_func(s0.transpose()*spins[m]);
          } else if (is_ghost) {
            dE = onsite_func(s0_new.transpose()*spins[j]) - onsite_func(s0.transpose()*spins[j]);
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
        set_spin(m, s_new);
      }
    }
  }
}

void Spin3DModel::metropolis_mutation() {
  nsteps++;
  acceptance = accepted/nsteps;
  if (acceptance > 0.5) {
    sigma = std::min(2., 1.01*sigma);
  } else {
    sigma = std::max(0.05, 0.99*sigma);
  }

  // Randomly generate mutation
  Eigen::Vector3d Gamma;
  Gamma << dist.sample(), dist.sample(), dist.sample();
  Eigen::Vector3d S2 = (spins[mut.i] + sigma*Gamma).normalized();


  // Store mutation for consideration
  mut.dS = S2 - spins[mut.i];
}

void Spin3DModel::generate_mutation() {
  if (cluster_update) {
    cluster_mutation();
  } else {
    mut.i = rand() % V;
    metropolis_mutation();
  }
}


void Spin3DModel::accept_mutation() {
  accepted++;
  return;
}

void Spin3DModel::reject_mutation() {
  if (!cluster_update) {
    set_spin(mut.i, spins[mut.i] - mut.dS);
  }
}

double Spin3DModel::onsite_energy(uint32_t i) const {
  if (cluster_update) {
    return onsite_func(s0.transpose()*spins[i]);
  } else {
    return onsite_func(spins[i]);
  }
}

double Spin3DModel::bond_energy(uint32_t i) const {
  double E = 0.;
  double dE;
  for (auto const &[j, b] : neighbors[i]) {
    if (b == GHOST) {
      continue;
    }
    dE = 0.5*bonds[b].bondfunc(spins[i], spins[j]);
    E += dE;
  }

  return E;
}

double Spin3DModel::energy() const {
  double E = 0;

  for (uint32_t i = 0; i < V; i++) {
    E += onsite_energy(i);
    E += bond_energy(i);
  }

  return E;
}

double Spin3DModel::energy_change() {
  if (cluster_update) {
    return -1.;
  }

  double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
  set_spin(mut.i, spins[mut.i] + mut.dS);
  double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

  return E2 - E1;
}

// Saves current spin configuration
void Spin3DModel::save_spins(const std::string& filename) {
  std::ofstream output_file;
  output_file.open(filename);
  output_file << N1 << "\t" << N2 << "\t" << N3 << "\t" << sl << "\n";
  Eigen::Vector3d S;
  for (uint32_t i = 0; i < V; i++) {
    if (cluster_update) {
      S = s0.transpose()*spins[i];
    } else {
      S = spins[i];
    }

    output_file << S[0] << "\t" << S[1] << "\t" << S[2];
    if (i < V-1) { 
      output_file << "\t";
    }
  }
  output_file.close();
}

void Spin3DModel::add_helicity_samples(dataframe::data_t &samples) const {
  std::vector<double> tterms = twist_derivatives();
  samples.emplace("d1E", tterms[0]);
  samples.emplace("d2E", tterms[1]);
  samples.emplace("d3E", tterms[2]);
  samples.emplace("d4E", tterms[3]);
}

void Spin3DModel::add_magnetization_samples(dataframe::data_t &samples) const {
  Eigen::Vector3d m = get_magnetization();
  samples.emplace("mx", m[0]);
  samples.emplace("my", m[1]);
  samples.emplace("mz", m[2]);
  samples.emplace("magnetization", m.norm());
}

dataframe::data_t Spin3DModel::take_samples() {
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

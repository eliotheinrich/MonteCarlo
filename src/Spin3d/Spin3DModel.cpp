#include "Spin3DModel.h"

#include <stack>
#include <stdexcept>

std::pair<std::vector<double>, std::vector<double>> split_complex_output(const std::vector<std::complex<double>>& complex_data) {
  size_t n = complex_data.size();
  std::vector<double> real(n);
  std::vector<double> imag(n);
  for (size_t i = 0; i < n; ++i) {
    real[i] = complex_data[i].real();
    imag[i] = complex_data[i].imag();
  }

  return std::make_pair(real, imag);
}

std::pair<std::vector<double>, std::vector<double>> fft2d_channel(std::vector<double>& x, std::vector<double>& y, std::vector<std::complex<double>>& input, int N, int s) {
  std::vector<std::complex<double>> output(N * N);
  finufft_opts opts;
  finufft_default_opts(&opts);
  opts.nthreads = 1;

  int ier = finufft2d1(s*N*N, &x[0], &y[0], &input[0], +1, 1e-6, N, N, &output[0], &opts);
  return split_complex_output(output);
}

std::pair<std::vector<double>, std::vector<double>> fft3d_channel(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<std::complex<double>>& input, int N, int s) {
  std::vector<std::complex<double>> output(N * N * N);
  finufft_opts opts;
  finufft_default_opts(&opts);
  opts.nthreads = 1;

  int ier = finufft3d1(s*N*N*N, &x[0], &y[0], &z[0], &input[0], +1, 1e-6, N, N, N, &output[0], &opts);
  return split_complex_output(output);
}

GaussianDist::GaussianDist(double mean, double std) {
  rd.seed(rand());
  gen = std::default_random_engine(rd());
  dist = std::normal_distribution<>(mean, std);
}

double GaussianDist::sample() {
  return dist(gen);
}

Spin3DModel::Spin3DModel(dataframe::ExperimentParams &params, uint32_t num_threads) : MonteCarloSimulator(params, num_threads), nsteps(0), accepted(0) {
  std::string mutation_str = dataframe::utils::get<std::string>(params, "mutation_type", "metropolis");
  if (mutation_str == "metropolis") {
    mutation_type = METROPOLIS;
  } else if (mutation_str == "adaptive_metropolis") {
    mutation_type = ADAPTIVE_METROPOLIS;
  } else if (mutation_str == "cluster") {
      mutation_type = CLUSTER;
  } else {
    throw std::runtime_error(fmt::format("Invalid mutation_type: {}", mutation_str));
  }

  if (mutation_type == METROPOLIS) {
    sigma = dataframe::utils::get<double>(params, "mutation_width", 2.0);
  } else if (mutation_type == ADAPTIVE_METROPOLIS) {
    sigma = 0.25;
  }

  acceptance = 0.5;

  sample_helicity      = dataframe::utils::get<int>(params, "sample_helicity",      false);
  sample_magnetization = dataframe::utils::get<int>(params, "sample_magnetization", true);
  sample_energy        = dataframe::utils::get<int>(params, "sample_energy",        true);
  sample_spins         = dataframe::utils::get<int>(params, "sample_spins",         false);

  dist = GaussianDist(0., 1.0);
}

void Spin3DModel::init(const Lattice<Spin3D>& lattice) {
  std::srand(randi());

  this->lattice = lattice;
  if (mutation_type == CLUSTER) {
    this->lattice.add_ghost();
  }

  V = lattice.system_size();

  s0 = Eigen::Matrix3d::Identity();

  randomize_spins();

  mut.i = randi(0, V);
}

void Spin3DModel::randomize_spins() {
  for (uint32_t i = 0; i < V; i++) {
    set_spin(i, Eigen::Vector3d::Random(3).normalized());
  }
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

  auto [R1s, R2s, R3s] = lattice.get_twist_matrices(alpha);

  for (auto const &[j, b] : lattice.neighbors[i]) {
    if (b == GHOST) {
      continue;
    }

    S2 = get_spin(j);

    E0  += lattice.bonds[b].bondfunc(S1, S2);

    E1  += lattice.bonds[b].bondfunc(S1, R1s[b]*S2);
    Em1 += lattice.bonds[b].bondfunc(S1, R1s[b].transpose()*S2);

    E2  += lattice.bonds[b].bondfunc(S1, R2s[b]*S2);
    Em2 += lattice.bonds[b].bondfunc(S1, R2s[b].transpose()*S2);

    E3  += lattice.bonds[b].bondfunc(S1, R3s[b]*S2);
    Em3 += lattice.bonds[b].bondfunc(S1, R3s[b].transpose()*S2);
  }

  double d1E = ( 1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha/2.;
  double d2E = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/pow(alpha, 2)/2.;
  double d3E = ( 1./8.*Em3 - 1.*Em2 + 13./8.*Em1 - 13./8.*E1 + 1.*E2 - 1./8.*E3)/pow(alpha, 3)/2.;
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
  Eigen::Vector3d M = Eigen::Vector3d::Zero();
  for (uint32_t i = 0; i < V; i++) {
    M += get_spin(i);
  }

  return M/V;
}

std::vector<double> Spin3DModel::correlation_function(uint32_t i, uint32_t a = 2, uint32_t b = 2) const {
  std::vector<double> Cij = std::vector<double>(V, 0.0); 

  auto spini = get_spin(i);
  for (uint32_t n1 = 0; n1 < lattice.dx.N; n1++) {
    for (uint32_t n2 = 0; n2 < lattice.dy.N; n2++) {
      for (uint32_t n3 = 0; n3 < lattice.dz.N; n3++) {
        for (uint32_t s = 0; s < lattice.num_sublattices; s++) {
          auto neighbor_opt = lattice.get_neighbor(i, n1, n2, n3, s);
          if (neighbor_opt) {
            uint32_t j = lattice.flat_idx(n1, n2, n3, s);
            Cij[j] = get_spin(neighbor_opt.value())[a] * spini[b];
          }
        }
      }
    }
  }
  return Cij;
}

std::vector<double> Spin3DModel::full_correlation_function(uint32_t i) const {
  std::vector<double> Cij = std::vector<double>(V, 0.0); 

  auto spini = get_spin(i);
  for (uint32_t n1 = 0; n1 < lattice.dx.N; n1++) {
    for (uint32_t n2 = 0; n2 < lattice.dy.N; n2++) {
      for (uint32_t n3 = 0; n3 < lattice.dz.N; n3++) {
        for (uint32_t s = 0; s < lattice.num_sublattices; s++) {
          auto neighbor_opt = lattice.get_neighbor(i, n1, n2, n3, s);
          if (neighbor_opt) {
            uint32_t j = lattice.flat_idx(n1, n2, n3, s);
            Cij[j] = get_spin(neighbor_opt.value()).dot(spini);
          }
        }
      }
    }
  }
  return Cij;
}

double Spin3DModel::skyrmion_density(uint32_t i) const {
  Eigen::Vector3d dSdX; dSdX << 0., 0., 0.;
  Eigen::Vector3d dSdY; dSdY << 0., 0., 0.;
  for (auto const &[j, b] : lattice.neighbors[i]) {
    if (b == GHOST) {
      continue;
    }

    auto bond = lattice.bonds[b];
    Eigen::Vector3d pos = lattice.position(bond.d1, bond.d2, bond.d3, bond.ds);
    dSdX += pos[0]*(lattice.spins[j] - lattice.spins[i]) + pos[1]*(lattice.spins[j] - lattice.spins[i]);
  }
  dSdX = dSdX/lattice.bonds.size();
  dSdY = dSdY/lattice.bonds.size();

  return lattice.spins[i].dot(dSdX.cross(dSdY));
}

std::vector<double> Spin3DModel::skyrmion_correlation_function(uint32_t i) const {
  std::vector<double> Cij = std::vector<double>(V, 0.0); 

  double Si = skyrmion_density(i);
  for (uint32_t n1 = 0; n1 < lattice.dx.N; n1++) {
    for (uint32_t n2 = 0; n2 < lattice.dy.N; n2++) {
      for (uint32_t n3 = 0; n3 < lattice.dz.N; n3++) {
        for (uint32_t s = 0; s < lattice.num_sublattices; s++) {
          auto neighbor_opt = lattice.get_neighbor(i, n1, n2, n3, s);
          if (neighbor_opt) {
            uint32_t j = lattice.flat_idx(n1, n2, n3, s);
            Cij[j] = skyrmion_density(neighbor_opt.value()) * Si;
          }
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
        s_new = R*lattice.spins[m];
      }

      for (auto const &[j, b] : lattice.neighbors[m]) {
        bool neighbor_is_ghost = (j == V);

        if (!s.count(j)) {
          double dE;

          if (neighbor_is_ghost) {
            dE = onsite_func(s0.transpose()*s_new) - onsite_func(s0.transpose()*lattice.spins[m]);
          } else if (is_ghost) {
            dE = onsite_func(s0_new.transpose()*lattice.spins[j]) - onsite_func(s0.transpose()*lattice.spins[j]);
          } else { // Normal bond
            dE = lattice.bonds[b].bondfunc(lattice.spins[j], s_new) - lattice.bonds[b].bondfunc(lattice.spins[j], lattice.spins[m]);
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
  // Basically random, with some small preference to nearby
  Eigen::Vector3d S2 = mutate_spin3d(lattice.spins[mut.i], dist, sigma);

  // Store mutation for consideration
  mut.dS = S2 - lattice.spins[mut.i];
}

void Spin3DModel::adaptive_metropolis_mutation() {
  nsteps++;
  acceptance = accepted/nsteps;
  if (acceptance > 0.5) {
    sigma = std::min(2., 1.01*sigma);
  } else {
    sigma = std::max(0.05, 0.99*sigma);
  }

  Eigen::Vector3d S2 = mutate_spin3d(lattice.spins[mut.i], dist, sigma);

  // Store mutation for consideration
  mut.dS = S2 - lattice.spins[mut.i];
}

void Spin3DModel::generate_mutation() {
  if (mutation_type == METROPOLIS) {
    mut.i = rand() % V;
    metropolis_mutation();
  } else if (mutation_type == ADAPTIVE_METROPOLIS) {
    mut.i = rand() % V;
    adaptive_metropolis_mutation();
  } else if (mutation_type == CLUSTER) {
    cluster_mutation();
  } 
}

void Spin3DModel::accept_mutation() {
  accepted++;
  return;
}

void Spin3DModel::reject_mutation() {
  if (mutation_type != CLUSTER) {
    set_spin(mut.i, lattice.spins[mut.i] - mut.dS);
  }
}

double Spin3DModel::onsite_energy(uint32_t i) const {
  if (mutation_type == CLUSTER) {
    return onsite_func(s0.transpose()*lattice.spins[i]);
  } else {
    return onsite_func(lattice.spins[i]);
  }
}

double Spin3DModel::bond_energy(uint32_t i) const {
  double E = 0.;
  for (auto const &[j, b] : lattice.neighbors[i]) {
    if (b == GHOST) {
      continue;
    }

    E += 0.5*lattice.bonds[b].bondfunc(lattice.spins[i], lattice.spins[j]);
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
  if (mutation_type == CLUSTER) {
    return -1.;
  }

  double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
  set_spin(mut.i, lattice.spins[mut.i] + mut.dS);
  double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

  return E2 - E1;
}

void Spin3DModel::add_magnetization_samples(dataframe::SampleMap &samples) const {
  Eigen::Vector3d m = get_magnetization();
  dataframe::utils::emplace(samples, "mx", m[0]);
  dataframe::utils::emplace(samples, "my", m[1]);
  dataframe::utils::emplace(samples, "mz", m[2]);

  dataframe::utils::emplace(samples, "mx_mag", std::abs(m[0]));
  dataframe::utils::emplace(samples, "my_mag", std::abs(m[1]));
  dataframe::utils::emplace(samples, "mz_mag", std::abs(m[2]));

  double mn = m.norm();
  dataframe::utils::emplace(samples, "magnetization", mn);
  dataframe::utils::emplace(samples, "magnetization_squared", mn*mn);
}

void Spin3DModel::add_helicity_samples(dataframe::SampleMap &samples) const {
  std::vector<double> tterms = twist_derivatives();
  dataframe::utils::emplace(samples, "d1E", tterms[0]);
  dataframe::utils::emplace(samples, "d2E", tterms[1]);
  dataframe::utils::emplace(samples, "d3E", tterms[2]);
  dataframe::utils::emplace(samples, "d4E", tterms[3]);
}

void Spin3DModel::add_spin_samples(dataframe::SampleMap& samples) const {
  std::vector<double> spins(3*V);
  for (size_t i = 0; i < V; i++) {
    auto spin = get_spin(i);
    spins[i]     = spin(0);
    spins[i+V]   = spin(1);
    spins[i+2*V] = spin(2);
  }

  dataframe::utils::emplace(samples, "spins", spins, {3, lattice.num_sublattices, lattice.dz.N, lattice.dy.N, lattice.dx.N});
}

dataframe::SampleMap Spin3DModel::take_samples() {
  dataframe::SampleMap samples;

  if (sample_energy) {
    double E = energy();
    dataframe::utils::emplace(samples, "energy", E);
    dataframe::utils::emplace(samples, "energy_squared", E*E);
  }

  if (sample_magnetization) {
    add_magnetization_samples(samples);
  }

  if (sample_helicity) {
    add_helicity_samples(samples);
  }

  if (sample_spins) {
    add_spin_samples(samples);
  }

  dataframe::utils::emplace(samples, "sigma", sigma);

  return samples;
}

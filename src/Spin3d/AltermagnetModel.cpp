#include "AltermagnetModel.h"

AltermagnetModel::AltermagnetModel(dataframe::ExperimentParams &params, uint32_t num_threads) : Spin3DModel(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");

  J1  = dataframe::utils::get<double>(params, "J1");
  J2  = dataframe::utils::get<double>(params, "J2");
  J2p = dataframe::utils::get<double>(params, "J2p");
  D1  = dataframe::utils::get<double>(params, "D1");
  D2  = dataframe::utils::get<double>(params, "D2");
  K   = dataframe::utils::get<double>(params, "K");
  B   = dataframe::utils::get<double>(params, "B", 0.0);
  T_i  = dataframe::utils::get<double>(params, "Ti", temperature);

  J1_i  = J1;
  J2_i  = J2;
  J2p_i = J2p;
  D1_i  = D1;
  D2_i  = D2;
  K_i   = K;
  B_i   = B;

  anneal = dataframe::utils::get<int>(params, "anneal", false);
  if (anneal) {
    T_f = dataframe::utils::get<double>(params, "Tf", T_i);
    J1_f  = dataframe::utils::get<double>(params, "J1_f", J1);
    J2_f  = dataframe::utils::get<double>(params, "J2_f", J2);
    J2p_f = dataframe::utils::get<double>(params, "J2p_f", J2p);
    D1_f  = dataframe::utils::get<double>(params, "D1_f", D1);
    D2_f  = dataframe::utils::get<double>(params, "D2_f", D2);
    K_f   = dataframe::utils::get<double>(params, "K_f", K);
    B_f   = dataframe::utils::get<double>(params, "B_f", 0.0);
  }

  initial_state = dataframe::utils::get<std::string>(params, "initial_state", "random");

  sample_sublattice_magnetization = dataframe::utils::get<int>(params, "sample_sublattice_magnetization", false);
  sample_structure_factor = dataframe::utils::get<int>(params, "sample_structure_factor", false);
  sample_staggered_structure_factor = dataframe::utils::get<int>(params, "sample_staggered_structure_factor", false);
}

void AltermagnetModel::init() {
  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfuncNN =
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return this->J1*S1.dot(S2);
    };

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfuncAx =
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return this->J2*S1.dot(S2);
    };

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfuncBx =
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return this->J2p*S1.dot(S2);
    };

  SiteFilter filterA = [](uint32_t i, uint32_t j, uint32_t k, uint32_t s) {
    return s == 0;
  };

  SiteFilter filterB = [](uint32_t i, uint32_t j, uint32_t k, uint32_t s) {
    return s == 1;
  };

  std::vector<SpinBond<Spin3D>> bonds = {
    SpinBond<Spin3D>( 1, 0, 0, 0, bondfuncAx, filterA), // n = 0
    SpinBond<Spin3D>(-1, 0, 0, 0, bondfuncAx, filterA), // n = 1
    SpinBond<Spin3D>( 0, 1, 0, 0, [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) { return this->J2p*S1.dot(S2) + this->D1*(S1[0] * S2[1] - S1[1] * S2[0]); }, filterA), // n = 2
    SpinBond<Spin3D>( 0,-1, 0, 0, [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) { return this->J2p*S1.dot(S2) - this->D1*(S1[0] * S2[1] - S1[1] * S2[0]); }, filterA), // n = 3

    SpinBond<Spin3D>( 1, 0, 0, 0, bondfuncBx, filterB), // n = 4
    SpinBond<Spin3D>(-1, 0, 0, 0, bondfuncBx, filterB), // n = 5
    SpinBond<Spin3D>( 0, 1, 0, 0, [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) { return this->J2 *S1.dot(S2) + this->D2*(S1[0] * S2[1] - S1[1] * S2[0]); }, filterB), // n = 6
    SpinBond<Spin3D>( 0,-1, 0, 0, [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) { return this->J2 *S1.dot(S2) - this->D2*(S1[0] * S2[1] - S1[1] * S2[0]); }, filterB), // n = 7

    SpinBond<Spin3D>( 0, 0, 0, 1, bondfuncNN),          // n = 8
    SpinBond<Spin3D>( 0,-1, 0, 1, bondfuncNN, filterA), // n = 9
    SpinBond<Spin3D>(-1,-1, 0, 1, bondfuncNN, filterA), // n = 10
    SpinBond<Spin3D>(-1, 0, 0, 1, bondfuncNN, filterA), // n = 11
    SpinBond<Spin3D>( 0, 1, 0, 1, bondfuncNN, filterB), // n = 12
    SpinBond<Spin3D>( 1, 1, 0, 1, bondfuncNN, filterB), // n = 13
    SpinBond<Spin3D>( 1, 0, 0, 1, bondfuncNN, filterB)  // n = 14
  };

  Eigen::Vector3d v1; v1 << 1., 0., 0.;
  Eigen::Vector3d v2; v2 << 0., 1., 0.;
  Eigen::Vector3d v3; v3 << 0., 0., 1.;
  Eigen::Vector3d vs; vs << 0.5, 0.5, 0.;

  LatticeDimension dx(N, BoundaryCondition::Periodic, v1);
  LatticeDimension dy(N, BoundaryCondition::Periodic, v2);
  LatticeDimension dz(1, BoundaryCondition::Periodic, v3);
  std::vector<Eigen::Vector3d> sublattice_vectors = {vs};
  Lattice<Spin3D> lattice_(dx, dy, dz, bonds, sublattice_vectors);
  Spin3DModel::init(lattice_);

  Spin3D z_polarized; z_polarized << 0.0, 0.0, 1.0;
  double sign = 1.0;
  if (randi() % 2) {
    sign = -1.0;
  }
  z_polarized *= sign;

  double noise = 0.2;

  if (initial_state == "g_afm") {
    for (size_t i = 0; i < V; i++) {
      auto idxs = lattice.tensor_idx(i);
      size_t s = idxs[3];
      if (s == 0) {
        lattice.spins[i] = mutate_spin3d( z_polarized, dist, noise);
      } else {
        lattice.spins[i] = mutate_spin3d(-z_polarized, dist, noise);
      }
    }
  } else if (initial_state == "fm") {
    for (size_t i = 0; i < V; i++) {
      lattice.spins[i] = mutate_spin3d(z_polarized, dist, noise);
    }
  } else if (initial_state == "stripe_afm_x") {
    for (size_t i = 0; i < V; i++) {
      auto idxs = lattice.tensor_idx(i);
      size_t s = idxs[3];
      double sublattice_sign = 1.0;
      if (s == 0) {
        sublattice_sign = -1.0;
      }

      size_t x = idxs[0];
      if (x % 2 == 0) {
        lattice.spins[i] = mutate_spin3d( sublattice_sign*z_polarized, dist, noise);
      } else {
        lattice.spins[i] = mutate_spin3d(-sublattice_sign*z_polarized, dist, noise);
      }
    }
  } else if (initial_state == "stripe_afm_y") {
    for (size_t i = 0; i < V; i++) {
      auto idxs = lattice.tensor_idx(i);
      size_t s = idxs[3];
      double sublattice_sign = 1.0;
      if (s == 0) {
        sublattice_sign = -1.0;
      }

      size_t y = idxs[0];
      if (y % 2 == 0) {
        lattice.spins[i] = mutate_spin3d( sublattice_sign*z_polarized, dist, noise);
      } else {
        lattice.spins[i] = mutate_spin3d(-sublattice_sign*z_polarized, dist, noise);
      }
    }
  }

}

// Does not include easy-axis anisotropy. Do not use!
void AltermagnetModel::over_relaxation_mutation() {
  Eigen::Vector3d H = Eigen::Vector3d::Zero();

  for (auto const &[j, n] : lattice.neighbors[mut.i]) {
    auto sj = get_spin(j);
    if (n == 0 || n == 1) {        // Ax
      H += J2*sj;         
    } else if (n == 2 || n == 3) { // Ay
      Eigen::Vector3d dv; dv << sj[1], -sj[0], 0;
      H += J2p*sj + D1*dv;
    } else if (n == 4 || n == 5) { // Bx
      H += J2p*sj;
    } else if (n == 6 || n == 7) { // By
      Eigen::Vector3d dv; dv << sj[1], -sj[0], 0;
      H += J2*sj + D2*dv;
    } else {                       // NN
      H += J1*sj;
    }
  }

  mut.dS = -2*get_spin(mut.i) + 2.*get_spin(mut.i).dot(H)/std::pow(H.norm(),2) * H;
}

double AltermagnetModel::onsite_func(const Eigen::Vector3d &S) const {
  double E = 0;
  E += K*S[2]*S[2];
  E += B*S[2];

  return E;
}

void AltermagnetModel::add_sublattice_magnetization_samples(dataframe::SampleMap &samples) const {
  Eigen::Vector3d M1 = Eigen::Vector3d::Zero();
  Eigen::Vector3d M2 = Eigen::Vector3d::Zero();

  for (size_t i = 0; i < V; i++) {
    auto idx = lattice.tensor_idx(i);
    if (idx[3] == 0) {
      M1 += get_spin(i);
    } else {
      M2 += get_spin(i);
    }
  }

  M1 = M1 / (V/2);
  M2 = M2 / (V/2);
  Eigen::Vector3d M = (M1 + M2)/2;
  Eigen::Vector3d neel = (M1 - M2)/2;

  double m1 = M1.norm();
  double m2 = M2.norm();
  double m  = M.norm();

  dataframe::utils::emplace(samples, "neelx", neel(0));
  dataframe::utils::emplace(samples, "neely", neel(1));
  dataframe::utils::emplace(samples, "neelz", neel(2));
  dataframe::utils::emplace(samples, "neel_mag", neel.norm());

  dataframe::utils::emplace(samples, "m1x", M1(0));
  dataframe::utils::emplace(samples, "m1y", M1(1));
  dataframe::utils::emplace(samples, "m1z", M1(2));

  dataframe::utils::emplace(samples, "m2x", M2(0));
  dataframe::utils::emplace(samples, "m2y", M2(1));
  dataframe::utils::emplace(samples, "m2z", M2(2));

  dataframe::utils::emplace(samples, "m1x_mag", std::abs(M1(0)));
  dataframe::utils::emplace(samples, "m1y_mag", std::abs(M1(1)));
  dataframe::utils::emplace(samples, "m1z_mag", std::abs(M1(2)));

  dataframe::utils::emplace(samples, "m2x_mag", std::abs(M2(0)));
  dataframe::utils::emplace(samples, "m2y_mag", std::abs(M2(1)));
  dataframe::utils::emplace(samples, "m2z_mag", std::abs(M2(2)));

  dataframe::utils::emplace(samples, "mx_mag", std::abs(M(0)));
  dataframe::utils::emplace(samples, "my_mag", std::abs(M(1)));
  dataframe::utils::emplace(samples, "mz_mag", std::abs(M(2)));

  dataframe::utils::emplace(samples, "magnetization", m);
  dataframe::utils::emplace(samples, "magnetization_squared", m*m);
  dataframe::utils::emplace(samples, "magnetization1", m1);
  dataframe::utils::emplace(samples, "magnetization2", m2);
}

void AltermagnetModel::add_structure_factor_samples(dataframe::SampleMap& samples, bool staggered=false) const {
  std::vector<std::complex<double>> Sx(V);
  std::vector<std::complex<double>> Sy(V);
  std::vector<std::complex<double>> Sz(V);

  std::vector<double> x(V);
  std::vector<double> y(V);

  for (size_t i = 0; i < V; i++) {
    auto pos = lattice.position(i);
    x[i] = 2.0 * M_PI * pos(0) / N - M_PI;
    y[i] = 2.0 * M_PI * pos(1) / N - M_PI;

    auto idxs = lattice.tensor_idx(i);
    auto spin = get_spin(i);
    if (staggered && idxs[3] == 0) {
      spin = -spin;
    }

    Sx[i] = spin[0];
    Sy[i] = spin[1];
    Sz[i] = spin[2];
  }

  auto to_mag = [](const std::vector<double>& real, const std::vector<double>& complex) {
    size_t n = real.size();
    std::vector<double> mag(n);
    for (size_t i = 0; i < n; i++) {
      mag[i] = real[i]*real[i] + complex[i]*complex[i];
    }
    return mag;
  };

  std::string prefix = staggered ? "s_" : "";

  auto [Sx1, Sx2] = fft2d_channel(x, y, Sx, N, 2);
  dataframe::utils::emplace(samples, fmt::format("{}Sx_real", prefix), Sx1,              {N, N});
  dataframe::utils::emplace(samples, fmt::format("{}Sx_imag", prefix), Sx2,              {N, N});
  dataframe::utils::emplace(samples, fmt::format("{}Sx",      prefix), to_mag(Sx1, Sx2), {N, N});

  auto [Sy1, Sy2] = fft2d_channel(x, y, Sy, N, 2);
  dataframe::utils::emplace(samples, fmt::format("{}Sy_real", prefix), Sy1,              {N, N});
  dataframe::utils::emplace(samples, fmt::format("{}Sy_imag", prefix), Sy2,              {N, N});
  dataframe::utils::emplace(samples, fmt::format("{}Sy", prefix),      to_mag(Sy1, Sy2), {N, N});

  auto [Sz1, Sz2] = fft2d_channel(x, y, Sz, N, 2);
  dataframe::utils::emplace(samples, fmt::format("{}Sz_real", prefix), Sz1,              {N, N});
  dataframe::utils::emplace(samples, fmt::format("{}Sz_imag", prefix), Sz2,              {N, N});
  dataframe::utils::emplace(samples, fmt::format("{}Sz", prefix),      to_mag(Sz1, Sz2), {N, N});
}

dataframe::SampleMap AltermagnetModel::take_samples() const {
  dataframe::SampleMap samples = Spin3DModel::take_samples();

  if (sample_sublattice_magnetization) {
    add_sublattice_magnetization_samples(samples);
  }

  if (sample_structure_factor) {
    add_structure_factor_samples(samples, false);
  }

  if (sample_staggered_structure_factor) {
    add_structure_factor_samples(samples, true);
  }

  return samples;
}


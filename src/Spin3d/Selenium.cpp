#include "Selenium.h"

SeleniumModel::SeleniumModel(dataframe::ExperimentParams &params, uint32_t num_threads) : Spin3DModel(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");

  J1  = dataframe::utils::get<double>(params, "J1");
  J2  = dataframe::utils::get<double>(params, "J2");
  J2p = dataframe::utils::get<double>(params, "J2p");
  K   = dataframe::utils::get<double>(params, "K");
  D1  = dataframe::utils::get<double>(params, "D1");
  D2  = dataframe::utils::get<double>(params, "D2");

  double D1l  = D1;
  double D2l  = D2;
  double J1l  = J1;
  double J2l  = J2;
  double J2pl = J2p;

  sample_sublattice_magnetization = dataframe::utils::get<int>(params, "sample_sublattice_magnetization", false);
  sample_structure_factor = dataframe::utils::get<int>(params, "sample_structure_factor", false);
  sample_staggered_structure_factor = dataframe::utils::get<int>(params, "sample_staggered_structure_factor", false);

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfuncNN =
    [J1l](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return J1l*S1.dot(S2);
    };

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfuncAx =
    [J2l](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return J2l*S1.dot(S2);
    };

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfuncBx =
    [J2pl](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return J2pl*S1.dot(S2);
    };

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfuncAy =
    [J2pl, D1l](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return J2pl*S1.dot(S2) + D1l*(S1[0] * S2[1] - S1[1] * S2[0]);
    };

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfuncBy =
    [J2l, D2l](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return J2l*S1.dot(S2) + D2l*(S1[0] * S2[1] - S1[1] * S2[0]);
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
    SpinBond<Spin3D>( 0, 1, 0, 0, bondfuncAy, filterA), // n = 2
    SpinBond<Spin3D>( 0,-1, 0, 0, bondfuncAy, filterA), // n = 3

    SpinBond<Spin3D>( 1, 0, 0, 0, bondfuncBx, filterB), // n = 4
    SpinBond<Spin3D>(-1, 0, 0, 0, bondfuncBx, filterB), // n = 5
    SpinBond<Spin3D>( 0, 1, 0, 0, bondfuncBy, filterB), // n = 6
    SpinBond<Spin3D>( 0,-1, 0, 0, bondfuncBy, filterB), // n = 7

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
  Lattice<Spin3D> lattice(dx, dy, dz, bonds, sublattice_vectors);
  Spin3DModel::init(lattice);
}

void SeleniumModel::generate_mutation() {
  if (cluster_update) {
    cluster_mutation(); 
  } else {
    mut.i = rand() % V;
    metropolis_mutation();
  }
}

// Does not include easy-axis anisotropy. Do not use!
void SeleniumModel::over_relaxation_mutation() {
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

double SeleniumModel::onsite_func(const Eigen::Vector3d &S) const {
  double E = 0;
  E += K*S[2]*S[2];

  //std::cout << "S = " << S << "\n";
  //std::cout << fmt::format("K = {}, onsite energy = {}\n", K, E);
  return E;
}

void SeleniumModel::add_sublattice_magnetization_samples(dataframe::SampleMap &samples) const {
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

void SeleniumModel::add_structure_factor_samples(dataframe::SampleMap& samples) const {
  std::vector<std::complex<double>> Sx(V);
  std::vector<std::complex<double>> Sy(V);
  std::vector<std::complex<double>> Sz(V);

  std::vector<double> x(V);
  std::vector<double> y(V);

  for (size_t i = 0; i < V; i++) {
    auto pos = lattice.position(i);
    x[i] = 2.0 * M_PI * pos(0) / N - M_PI;
    y[i] = 2.0 * M_PI * pos(1) / N - M_PI;

    auto spin = get_spin(i);
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

  auto [Sx1, Sx2] = fft2d_channel(x, y, Sx, N, 2);
  dataframe::utils::emplace(samples, "Sx_real", Sx1, {N, N});
  dataframe::utils::emplace(samples, "Sx_imag", Sx2, {N, N});
  dataframe::utils::emplace(samples, "Sx", to_mag(Sx1, Sx2), {N, N});

  auto [Sy1, Sy2] = fft2d_channel(x, y, Sy, N, 2);
  dataframe::utils::emplace(samples, "Sy_real", Sy1, {N, N});
  dataframe::utils::emplace(samples, "Sy_imag", Sy2, {N, N});
  dataframe::utils::emplace(samples, "Sy", to_mag(Sy1, Sy2), {N, N});

  auto [Sz1, Sz2] = fft2d_channel(x, y, Sz, N, 2);
  dataframe::utils::emplace(samples, "Sz_real", Sz1, {N, N});
  dataframe::utils::emplace(samples, "Sz_imag", Sz2, {N, N});
  dataframe::utils::emplace(samples, "Sz", to_mag(Sz1, Sz2), {N, N});
}

void SeleniumModel::add_staggered_structure_factor_samples(dataframe::SampleMap& samples) const {
  std::vector<std::complex<double>> Sx(V);
  std::vector<std::complex<double>> Sy(V);
  std::vector<std::complex<double>> Sz(V);

  std::vector<double> x(V);
  std::vector<double> y(V);

  for (size_t i = 0; i < V; i++) {
    auto pos = lattice.position(i);
    x[i] = 2.0 * M_PI * pos(0) / N - M_PI;
    y[i] = 2.0 * M_PI * pos(1) / N - M_PI;

    auto spin = get_spin(i);

    auto idxs = lattice.tensor_idx(i);
    if (idxs[3] == 1) {
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

  auto [Sx1, Sx2] = fft2d_channel(x, y, Sx, N, 2);
  dataframe::utils::emplace(samples, "_Sx_real", Sx1, {N, N});
  dataframe::utils::emplace(samples, "_Sx_imag", Sx2, {N, N});
  dataframe::utils::emplace(samples, "_Sx", to_mag(Sx1, Sx2), {N, N});

  auto [Sy1, Sy2] = fft2d_channel(x, y, Sy, N, 2);
  dataframe::utils::emplace(samples, "_Sy_real", Sy1, {N, N});
  dataframe::utils::emplace(samples, "_Sy_imag", Sy2, {N, N});
  dataframe::utils::emplace(samples, "_Sy", to_mag(Sy1, Sy2), {N, N});

  auto [Sz1, Sz2] = fft2d_channel(x, y, Sz, N, 2);
  dataframe::utils::emplace(samples, "_Sz_real", Sz1, {N, N});
  dataframe::utils::emplace(samples, "_Sz_imag", Sz2, {N, N});
  dataframe::utils::emplace(samples, "_Sz", to_mag(Sz1, Sz2), {N, N});
}

dataframe::SampleMap SeleniumModel::take_samples() {
  dataframe::SampleMap samples = Spin3DModel::take_samples();

  if (sample_sublattice_magnetization) {
    add_sublattice_magnetization_samples(samples);
  }

  if (sample_structure_factor) {
    add_structure_factor_samples(samples);
  }

  if (sample_staggered_structure_factor) {
    add_staggered_structure_factor_samples(samples);
  }

  return samples;
}


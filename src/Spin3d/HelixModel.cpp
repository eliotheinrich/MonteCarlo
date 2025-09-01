#include "HelixModel.h"

HelixModel::HelixModel(dataframe::ExperimentParams &params, uint32_t num_threads) : Spin3DModel(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");

  J1  = dataframe::utils::get<double>(params, "J1");
  J2  = dataframe::utils::get<double>(params, "J2");
  D1  = dataframe::utils::get<double>(params, "D1");
  D2  = dataframe::utils::get<double>(params, "D2");

  anneal = dataframe::utils::get<int>(params, "anneal", 0);
  min_temp = dataframe::utils::get<double>(params, "min_temp", temperature);
  max_temp = dataframe::utils::get<double>(params, "max_temp", temperature);

  sample_structure_factor = dataframe::utils::get<int>(params, "sample_structure_factor", false);

  Eigen::Vector3d xhat; xhat << 1.0, 0.0, 0.0;
  Eigen::Vector3d yhat; yhat << 0.0, 1.0, 0.0;
  Eigen::Vector3d zhat; zhat << 0.0, 0.0, 1.0;

  auto bondfunc = [](double J, double D, const Eigen::Vector3d &Dvec) {
    return [J, D, Dvec](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return -J * S1.dot(S2) - D * Dvec.dot(S1.cross(S2));
    };
  };

  std::vector<SpinBond<Spin3D>> bonds = {
    SpinBond<Spin3D>( 1, 0, 0, 0, bondfunc(J1, D1, xhat)),  // n = 0
    SpinBond<Spin3D>(-1, 0, 0, 0, bondfunc(J1,-D1, xhat)),  // n = 1
    SpinBond<Spin3D>( 0, 1, 0, 0, bondfunc(J1, D1, yhat)),  // n = 2
    SpinBond<Spin3D>( 0,-1, 0, 0, bondfunc(J1,-D1, yhat)),  // n = 3
    SpinBond<Spin3D>( 0, 0, 1, 0, bondfunc(J1, D1, zhat)),  // n = 4
    SpinBond<Spin3D>( 0, 0,-1, 0, bondfunc(J1,-D1, zhat)),  // n = 5
                                              
    SpinBond<Spin3D>( 2, 0, 0, 0, bondfunc(J2, D2, xhat)),  // n = 6
    SpinBond<Spin3D>(-2, 0, 0, 0, bondfunc(J2,-D2, xhat)),  // n = 7
    SpinBond<Spin3D>( 0, 2, 0, 0, bondfunc(J2, D2, yhat)),  // n = 8
    SpinBond<Spin3D>( 0,-2, 0, 0, bondfunc(J2,-D2, yhat)),  // n = 9
    SpinBond<Spin3D>( 0, 0, 2, 0, bondfunc(J2, D2, zhat)),  // n = 10
    SpinBond<Spin3D>( 0, 0,-2, 0, bondfunc(J2,-D2, zhat)),  // n = 11
  };

  Eigen::Vector3d v1; v1 << 1., 0., 0.;
  Eigen::Vector3d v2; v2 << 0., 1., 0.;
  Eigen::Vector3d v3; v3 << 0., 0., 1.;

  LatticeDimension dx(N, BoundaryCondition::Periodic, v1);
  LatticeDimension dy(N, BoundaryCondition::Periodic, v2);
  LatticeDimension dz(N, BoundaryCondition::Periodic, v3);
  Lattice<Spin3D> lattice(dx, dy, dz, bonds);
  Spin3DModel::init(lattice);
}

void HelixModel::add_structure_factor_samples(dataframe::SampleMap& samples) const {
  std::vector<std::complex<double>> Sx(V);
  std::vector<std::complex<double>> Sy(V);
  std::vector<std::complex<double>> Sz(V);

  std::vector<double> x(V);
  std::vector<double> y(V);
  std::vector<double> z(V);

  for (size_t i = 0; i < V; i++) {
    auto pos = lattice.position(i);
    x[i] = 2.0 * M_PI * pos(0) / N - M_PI;
    y[i] = 2.0 * M_PI * pos(1) / N - M_PI;
    z[i] = 2.0 * M_PI * pos(2) / N - M_PI;

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

  auto [Sx1, Sx2] = fft3d_channel(x, y, z, Sx, N, 1);
  dataframe::utils::emplace(samples, "Sx_real", Sx1, {N, N, N});
  dataframe::utils::emplace(samples, "Sx_imag", Sx2, {N, N, N});
  dataframe::utils::emplace(samples, "Sx", to_mag(Sx1, Sx2), {N, N, N});

  auto [Sy1, Sy2] = fft3d_channel(x, y, z, Sy, N, 1);
  dataframe::utils::emplace(samples, "Sy_real", Sy1, {N, N, N});
  dataframe::utils::emplace(samples, "Sy_imag", Sy2, {N, N, N});
  dataframe::utils::emplace(samples, "Sy", to_mag(Sy1, Sy2), {N, N, N});

  auto [Sz1, Sz2] = fft3d_channel(x, y, z, Sz, N, 1);
  dataframe::utils::emplace(samples, "Sz_real", Sz1, {N, N, N});
  dataframe::utils::emplace(samples, "Sz_imag", Sz2, {N, N, N});
  dataframe::utils::emplace(samples, "Sz", to_mag(Sz1, Sz2), {N, N, N});
}

dataframe::SampleMap HelixModel::take_samples() const {
  dataframe::SampleMap samples = Spin3DModel::take_samples();

  if (sample_structure_factor) {
    add_structure_factor_samples(samples);
  }

  return samples;
}


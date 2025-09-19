#include "AltermagnetModel.h"

AltermagnetModel::AltermagnetModel(Params &params, uint32_t num_threads) : Spin3DModel(params, num_threads) {
  N = get<int>(params, "system_size");

  J1  = get<double>(params, "J1");
  J2  = get<double>(params, "J2");
  J2p = get<double>(params, "J2p");
  D1  = get<double>(params, "D1");
  D2  = get<double>(params, "D2");
  K   = get<double>(params, "K");
  B   = get<double>(params, "B", 0.0);

  J1_i  = J1;
  J2_i  = J2;
  J2p_i = J2p;
  D1_i  = D1;
  D2_i  = D2;
  K_i   = K;
  B_i   = B;
  T_i = temperature;

  anneal = get<int>(params, "anneal", false);
  if (anneal) {
    T_f = get<double>(params, "T_f", T_i);
    J1_f  = get<double>(params, "J1_f", J1);
    J2_f  = get<double>(params, "J2_f", J2);
    J2p_f = get<double>(params, "J2p_f", J2p);
    D1_f  = get<double>(params, "D1_f", D1);
    D2_f  = get<double>(params, "D2_f", D2);
    K_f   = get<double>(params, "K_f", K);
    B_f   = get<double>(params, "B_f", 0.0);
  }

  initial_state = get<std::string>(params, "initial_state", "random");
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

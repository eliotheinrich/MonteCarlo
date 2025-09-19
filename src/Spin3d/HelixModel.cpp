#include "HelixModel.h"

HelixModel::HelixModel(Params &params, uint32_t num_threads) : Spin3DModel(params, num_threads) {
  N = get<int>(params, "system_size");

  J1  = get<double>(params, "J1");
  J2  = get<double>(params, "J2");
  D1  = get<double>(params, "D1");
  D2  = get<double>(params, "D2");

  anneal = get<int>(params, "anneal", 0);
  min_temp = get<double>(params, "min_temp", temperature);
  max_temp = get<double>(params, "max_temp", temperature);

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

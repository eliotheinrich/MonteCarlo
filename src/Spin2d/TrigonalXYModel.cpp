#include "TrigonalXYModel.h"
#include <functional>

TrigonalXYModel::TrigonalXYModel(dataframe::ExperimentParams &params, uint32_t num_threads) : Spin2DModel(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");
  L = dataframe::utils::get<int>(params, "layers", DEFAULT_LAYERS);
  J = dataframe::utils::get<double>(params, "J");
  A = dataframe::utils::get<double>(params, "A");

  double Jt = J;
  std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)> dotfunc = 
    [Jt](const Eigen::Vector2d &S1, const Eigen::Vector2d &S2) {
      return -Jt*S1.dot(S2);
    };

  std::vector<SpinBond<Spin2D>> bonds = {
    SpinBond<Spin2D>( 1, 0, 0, 0, dotfunc),
    SpinBond<Spin2D>(-1, 0, 0, 0, dotfunc),
    SpinBond<Spin2D>( 0, 1, 0, 0, dotfunc),
    SpinBond<Spin2D>( 0,-1, 0, 0, dotfunc),
    SpinBond<Spin2D>( 1,-1, 0, 0, dotfunc),
    SpinBond<Spin2D>(-1, 1, 0, 0, dotfunc)
  };

  Eigen::Vector3d v1; v1 << 1., 0., 0.;
  Eigen::Vector3d v2; v2 << 0.5, std::sqrt(3)/2., 0.;
  Eigen::Vector3d v3; v3 << 0., 0., 1.;

  LatticeDimension dx(N, BoundaryCondition::Periodic, v1);
  LatticeDimension dy(N, BoundaryCondition::Periodic, v2);
  LatticeDimension dz(L, BoundaryCondition::Periodic, v3);
  
  Lattice<Spin2D> lattice(dx, dy, dz, bonds);
  Spin2DModel::init(lattice);

  mut_mode = 0;
}


std::vector<double> TrigonalXYModel::vorticity() const {
  double v1 = 0;
  double v2 = 0;

  std::vector<std::vector<std::vector<double>>> phi = std::vector<std::vector<std::vector<double>>>(N,
                                                                  std::vector<std::vector<double>>(N,
                                                                              std::vector<double>(L)));

  for (uint32_t n1 = 0; n1 < N; n1++) {
    for (uint32_t n2 = 0; n2 < N; n2++) {
      for (uint32_t n3 = 0; n3 < L; n3++) {
        uint32_t i = lattice.flat_idx(n1, n2, n3, 0);
        phi[n1][n2][n3] = atan2(get_spin(i)[1], get_spin(i)[0]);
      }
    }
  }

  for (uint32_t n1 = 0; n1 < N; n1++) {
    for (uint32_t n2 = 0; n2 < N; n2++) {
      for (uint32_t n3 = 0; n3 < L; n3++) {
        double p1 = phi[n1][n2][n3]; 
        double p2 = phi[(n1+1)%N][n2][n3];
        double p3 = phi[(n1+1)%N][(n2+1)%N][n3];
        double w = std::arg(exp(std::complex<double>(0., p2 - p1))) + std::arg(exp(std::complex<double>(0., p3 - p2)))
                 + std::arg(exp(std::complex<double>(0., p1 - p3)));
        if (w > 0) { 
          v1 += w; 
        } else { 
          v2 += w; 
        }
      }
    }
  }

  return std::vector<double>{v1/(2*M_PI*N*N*L), v2/(2*M_PI*N*N*L)};
}

double TrigonalXYModel::onsite_func(const Eigen::Vector2d &S) const {
  double phi = atan2(S[1], S[0]);
  return A*std::cos(6*phi);
}

void TrigonalXYModel::over_relaxation_mutation() {
  Eigen::Vector2d H; H << 0., 0.;
  for (uint32_t n = 0; n < lattice.bonds.size(); n++) {
    uint32_t j = lattice.neighbors[mut.i][n].first;
    H -= J*get_spin(j);
  }

  mut.dS = -2*get_spin(mut.i) + 2.*get_spin(mut.i).dot(H)/std::pow(H.norm(),2) * H;
}

void TrigonalXYModel::generate_mutation() {
  if (cluster_update) {
    cluster_mutation(); 
  } else {
    mut.i = rand() % V;

    if (mut.i == 0) {
      mut_mode++;
    }

    if (mut_mode < 10) {
      over_relaxation_mutation();
    } else if (mut_mode < 14) {
      metropolis_mutation();
    } else {
      metropolis_mutation();
      mut_mode = 0;
    }
  }
}

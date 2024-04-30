#include "XXZHeis.h"
#include <functional>


XXZHeis::XXZHeis(dataframe::Params &params, uint32_t num_threads) : Spin3DModel(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");
  L = dataframe::utils::get<int>(params, "layers", DEFAULT_LAYERS);

  J = dataframe::utils::get<double>(params, "J");
  K = dataframe::utils::get<double>(params, "K");
  A = dataframe::utils::get<double>(params, "A");

  double Kt = K;
  double Jt = J;
  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc = 
    [Jt, Kt](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return -Jt*S1.dot(S2) + Kt*S1[2]*S2[2];
    };



  Eigen::Vector3d v1; v1 << 1.,0.,0.;
  Eigen::Vector3d v2; v2 << 0.,1.,0.;
  add_bond(1,0,0,0,   v1, bondfunc);
  add_bond(-1,0,0,0, -v1, bondfunc);
  add_bond(0,1,0,0,   v2, bondfunc);
  add_bond(0,-1,0,0, -v2, bondfunc);

  Spin3DModel::init(1, N, N, L);
}

std::vector<double> XXZHeis::vorticity() const {
  double v1 = 0;
  double v2 = 0;

  std::vector<std::vector<std::vector<double>>> phi = std::vector<std::vector<std::vector<double>>>(N,
                                                                  std::vector<std::vector<double>>(N,
                                                                              std::vector<double>(L)));

  int i;
  for (uint32_t n1 = 0; n1 < N; n1++) {
    for (uint32_t n2 = 0; n2 < N; n2++) {
      for (uint32_t n3 = 0; n3 < L; n3++) {
        i = flat_idx(n1, n2, n3, 0);
        phi[n1][n2][n3] = atan2(get_spin(i)[1], get_spin(i)[0]);
      }
    }
  }

  double p1; double p2; double p3; double p4;
  double w;
  for (uint32_t n1 = 0; n1 < N; n1++) {
    for (uint32_t n2 = 0; n2 < N; n2++) {
      for (uint32_t n3 = 0; n3 < L; n3++) {
        p1 = phi[n1][n2][n3]; p2 = phi[(n1+1)%N][n2][n3];
        p3 = phi[(n1+1)%N][(n2+1)%N][n3]; p4 = phi[n1][(n2+1)%N][n3];
        w = arg(exp(std::complex<double>(0., p2 - p1))) + arg(exp(std::complex<double>(0., p3 - p2)))
          + arg(exp(std::complex<double>(0., p4 - p3))) + arg(exp(std::complex<double>(0., p1 - p4)));
        if (w > 0) { 
          v1 += w; 
        } else { 
          v2 += w; 
        }
      }
    }
  }

  return std::vector<double>{v1/(2*PI*N*N*L), v2/(2*PI*N*N*L)};
}

double XXZHeis::onsite_func(const Eigen::Vector3d &S) const {
  // Onsite interactions
  return A*S[2]*S[2];
}

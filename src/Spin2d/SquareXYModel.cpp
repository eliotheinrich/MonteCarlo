#include "SquareXYModel.h"
#include <complex>
#include <functional>

SquareXYModel::SquareXYModel(dataframe::ExperimentParams &params, uint32_t num_threads) : Spin2DModel(params, num_threads) {
    N = dataframe::utils::get<int>(params, "system_size");
    L = dataframe::utils::get<int>(params, "layers", DEFAULT_LAYERS);

    J = dataframe::utils::get<double>(params, "J");
    B = dataframe::utils::get<double>(params, "B");
    Bp = dataframe::utils::get<double>(params, "Bp");
    Bx = B*cos(Bp);
    By = B*sin(Bp);

    double Jt = J;
    std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)> bondfunc = 
    [Jt](const Eigen::Vector2d &S1, const Eigen::Vector2d &S2) {
        return -Jt*S1.dot(S2);
    };

    if (!cluster_update) {
      mut_mode = 0;
    }

    std::vector<SpinBond<Spin2D>> bonds = {
      SpinBond<Spin2D>( 1, 0, 0, 0, bondfunc),
      SpinBond<Spin2D>(-1, 0, 0, 0, bondfunc),
      SpinBond<Spin2D>( 0, 1, 0, 0, bondfunc),
      SpinBond<Spin2D>( 0,-1, 0, 0, bondfunc),
    };

    Eigen::Vector3d v1; v1 << 1.,0.,0.;
    Eigen::Vector3d v2; v2 << 0.,1.,0.;
    Eigen::Vector3d v3; v3 << 0.,0.,1.;

    LatticeDimension dx(N, BoundaryCondition::Periodic, v1);
    LatticeDimension dy(N, BoundaryCondition::Periodic, v2);
    LatticeDimension dz(L, BoundaryCondition::Periodic, v3);

    Lattice<Spin2D> lattice(dx, dy, dz, bonds);
    Spin2DModel::init(lattice);
}

std::vector<std::pair<size_t, bool>> SquareXYModel::get_vortices() const {
  std::vector<std::vector<std::vector<double>>> phi = std::vector<std::vector<std::vector<double>>>(N,
      std::vector<std::vector<double>>(N,
        std::vector<double>(L)));

  for (uint32_t n1 = 0; n1 < N; n1++) {
    for (uint32_t n2 = 0; n2 < N; n2++) {
      for (uint32_t n3 = 0; n3 < L; n3++) {
        uint32_t i = lattice.flat_idx(n1, n2, n3, 0);
        phi[n1][n2][n3] = std::atan2(get_spin(i)[1], get_spin(i)[0]);
      }
    }
  }

  std::vector<std::pair<size_t, bool>> vortex_indices;
  for (uint32_t n1 = 0; n1 < N; n1++) {
    for (uint32_t n2 = 0; n2 < N; n2++) {
      for (uint32_t n3 = 0; n3 < L; n3++) {
        size_t i = lattice.flat_idx(n1, n2, n3, 0);
        double p1 = phi[n1][n2][n3]; 
        double p2 = phi[(n1+1)%N][n2][n3];
        double p3 = phi[(n1+1)%N][(n2+1)%N][n3]; 
        double p4 = phi[n1][(n2+1)%N][n3];
        double w = (arg(exp(std::complex<double>(0., p2 - p1))) + arg(exp(std::complex<double>(0., p3 - p2)))
                  + arg(exp(std::complex<double>(0., p4 - p3))) + arg(exp(std::complex<double>(0., p1 - p4))))/(2*M_PI);
        if (w > 1e-4) { 
          vortex_indices.push_back(std::make_pair(i, true));
        } else if (w < -1e-4) { 
          vortex_indices.push_back(std::make_pair(i, false));
        }
      }
    }
  }

  return vortex_indices; 
}

std::vector<double> SquareXYModel::vorticity() const {
  auto vortex_locations = get_vortices();

  double v1 = 0.0;
  double v2 = 0.0;
  for (const auto& [idx, v] : vortex_locations) {
    if (v) {
      v1 += 1;
    } else {
      v2 += 1;
    }
  }

  return std::vector<double>{v1/(N*N*L), v2/(N*N*L)};
}

double SquareXYModel::p(uint32_t i) const {
  return std::atan2(get_spin(i)[1], get_spin(i)[0]);
}

double SquareXYModel::e1() const {
  double s = 0;
  for (uint32_t i = 0; i < V; i++) {
    s += std::cos(p(i) - p(lattice.neighbors[i][0].first));
  }
  return s/V;
}

double SquareXYModel::e2() const {
  double s = 0;
  for (uint32_t i = 0; i < V; i++) {
    s += std::sin(p(i) - p(lattice.neighbors[i][0].first));
  }
  return s/V;
}

double SquareXYModel::U2() const {
    return e1() - V/temperature*std::pow(e2(), 2);
}

std::vector<double> SquareXYModel::twist_stiffness() const {
  // Returns the first and second derivative in response to a phase twist
  double E0 = 0.;
  double E1 = 0.;
  double E2 = 0.;
  double E3 = 0.;
  double Em1 = 0.;
  double Em2 = 0.;
  double Em3 = 0.;

  auto [R1s, R2s, R3s] = lattice.get_twist_matrices(alpha);

  for (uint32_t i = 0; i < V; i++) {
    int j = lattice.neighbors[i][0].first;

    Spin2D S1 = get_spin(i);
    Spin2D S2 = get_spin(j);

    E0  += lattice.bonds[0].bondfunc(S1, S2);

    E1  += lattice.bonds[0].bondfunc(S1, R1s[0]*S2);
    Em1 += lattice.bonds[0].bondfunc(S1, R1s[0].transpose()*S2);

    E2  += lattice.bonds[0].bondfunc(S1, R2s[0]*S2);
    Em2 += lattice.bonds[0].bondfunc(S1, R2s[0].transpose()*S2);

    E3  += lattice.bonds[0].bondfunc(S1, R3s[0]*S2);
    Em3 += lattice.bonds[0].bondfunc(S1, R3s[0].transpose()*S2);
  }

  // Compute derivates from finite difference
  double d1E = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha;
  double d2E = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/pow(alpha, 2);
  double d3E = (1./8.*Em3 - 1.*Em2 + 13./8.*Em1 - 13./8.*E1 + 1.*E2 - 1./8.*E3)/pow(alpha, 3);
  double d4E = (-1./6.*Em3 + 2.*Em2 - 13./2.*Em1 + 28./3.*E0 - 13./2.*E1 + 2.*E2 - 1./6.*E3)/pow(alpha, 4);

  return std::vector<double>{d4E, d3E, d1E, pow(d2E,2), d2E, d1E*d3E, pow(d1E,2)*d2E, 
    d1E*d2E, pow(d1E,2), pow(d1E,4), pow(d1E,3), e1(), pow(e2(),4), U2()};
}

double SquareXYModel::onsite_func(const Spin2D &S) const {
  // Onsite interactions
  return -Bx*S[0] - By*S[1];
}

void SquareXYModel::over_relaxation_mutation() {
  Eigen::Vector2d H; H << 0., 0.;
  for (uint32_t n = 0; n < lattice.bonds.size(); n++) {
    int j = lattice.neighbors[mut.i][n].first;
    H -= J*get_spin(j);
  }

  mut.dS = -2*get_spin(mut.i) + 2.*get_spin(mut.i).dot(H)/std::pow(H.norm(),2) * H;
}

void SquareXYModel::generate_mutation() {
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

Texture SquareXYModel::get_texture() const {
  Texture texture(N, N);
  for (size_t i = 0; i < V; i++) {
    auto idxs = lattice.tensor_idx(i);
    size_t x = idxs[0];
    size_t y = idxs[1];
    double phi = std::atan2(get_spin(i)[1], get_spin(i)[0]);
    float c = (std::sin(phi) + 1.0)/2.0;
    texture.set(x, y, {c, c, c, 1.0});
  }

  auto vortices = get_vortices();
  for (const auto& [i, v] : vortices) {
    auto idxs = lattice.tensor_idx(i);
    size_t x = idxs[0];
    size_t y = idxs[1];

    if (v) {
      texture.set(x, y, {1.0, 0.0, 0.0, 1.0});
      texture.set(mod(x+1, N), y, {1.0, 0.0, 0.0, 1.0});
      texture.set(mod(x-1, N), y, {1.0, 0.0, 0.0, 1.0});
      texture.set(x, mod(y+1, N), {1.0, 0.0, 0.0, 1.0});
      texture.set(x, mod(y-1, N), {1.0, 0.0, 0.0, 1.0});
    } else {
      texture.set(x, y, {0.0, 0.0, 1.0, 1.0});
      texture.set(mod(x+1, N), y, {0.0, 0.0, 1.0, 1.0});
      texture.set(mod(x-1, N), y, {0.0, 0.0, 1.0, 1.0});
    }
  }

  return texture;
}

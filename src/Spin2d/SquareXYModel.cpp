#include "SquareXYModel.h"
#include <complex>
#include <functional>

SquareXYModel::SquareXYModel(dataframe::Params &params, uint32_t num_threads) : Spin2DModel(params, num_threads) {
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

    Eigen::Vector3d v1; v1 << 1.,0.,0.;
    Eigen::Vector3d v2; v2 << 0.,1.,0.;
    add_bond(1,0,0,0,   v1, bondfunc);
    add_bond(-1,0,0,0, -v1, bondfunc);
    add_bond(0,1,0,0,   v2, bondfunc);
    add_bond(0,-1,0,0, -v2, bondfunc);

    Spin2DModel::init(1, N, N, L);
}

std::vector<double> SquareXYModel::vorticity() const {
    double v1 = 0;
    double v2 = 0;

    std::vector<std::vector<std::vector<double>>> phi = std::vector<std::vector<std::vector<double>>>(N,
                                                                    std::vector<std::vector<double>>(N,
                                                                                std::vector<double>(L)));

    uint32_t i;
    for (uint32_t n1 = 0; n1 < N; n1++) {
        for (uint32_t n2 = 0; n2 < N; n2++) {
            for (uint32_t n3 = 0; n3 < L; n3++) {
                i = flat_idx(n1, n2, n3, 0);
                phi[n1][n2][n3] = 0.;
                phi[n1][n2][n3] = std::atan2(get_spin(i)[1], get_spin(i)[0]);
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

double SquareXYModel::p(uint32_t i) const {
    return std::atan2(get_spin(i)[1], get_spin(i)[0]);
}

double SquareXYModel::e1() const {
    double s = 0;
    for (uint32_t i = 0; i < V; i++) {
        s += std::cos(p(i) - p(neighbors[i][0].first));
    }
    return s/V;
}

double SquareXYModel::e2() const {
    double s = 0;
    for (uint32_t i = 0; i < V; i++) {
        s += std::sin(p(i) - p(neighbors[i][0].first));
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

    for (uint32_t i = 0; i < V; i++) {
        int j = neighbors[i][0].first;

        Eigen::Vector2d S1 = get_spin(i);
        Eigen::Vector2d S2 = get_spin(j);

        E0 += bonds[0].bondfunc(S1, S2);

        E1 += bonds[0].bondfunc(S1, R1s[0]*S2);
        Em1 += bonds[0].bondfunc(S1, R1s[0].transpose()*S2);

        E2 += bonds[0].bondfunc(S1, R2s[0]*S2);
        Em2 += bonds[0].bondfunc(S1, R2s[0].transpose()*S2);

        E3 += bonds[0].bondfunc(S1, R3s[0]*S2);
        Em3 += bonds[0].bondfunc(S1, R3s[0].transpose()*S2);
    }

    // Compute derivates from finite difference
    double d1E = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha;
    double d2E = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/pow(alpha, 2);
    double d3E = (1./8.*Em3 - 1.*Em2 + 13./8.*Em1 - 13./8.*E1 + 1.*E2 - 1./8.*E3)/pow(alpha, 3);
    double d4E = (-1./6.*Em3 + 2.*Em2 - 13./2.*Em1 + 28./3.*E0 - 13./2.*E1 + 2.*E2 - 1./6.*E3)/pow(alpha, 4);
    
    return std::vector<double>{d4E, d3E, d1E, pow(d2E,2), d2E, d1E*d3E, pow(d1E,2)*d2E, 
                               d1E*d2E, pow(d1E,2), pow(d1E,4), pow(d1E,3), e1(), pow(e2(),4), U2()};
}

double SquareXYModel::onsite_func(const Eigen::Vector2d &S) const {
    // Onsite interactions
    return -Bx*S[0] - By*S[1];
}

void SquareXYModel::over_relaxation_mutation() {
    Eigen::Vector2d H; H << 0., 0.;
    for (uint32_t n = 0; n < bonds.size(); n++) {
        int j = neighbors[mut.i][n].first;
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

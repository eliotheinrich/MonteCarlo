#include "SquareXYModel.h"
#include <iostream>
#include <complex>
#include <functional>
#include <string>

SquareXYModel::SquareXYModel(Params &params) : Spin2DModel(params) {
    this->N = params.get<int>("system_size");
    this->L = params.get<int>("layers", DEFAULT_LAYERS);

    Spin2DModel::init_params(1, N, N, L);

    this->J = params.get<float>("J");
    this->B = params.get<float>("B");
    this->Bp = params.get<float>("Bp");
    this->Bx = this->B*cos(this->Bp);
    this->By = this->B*sin(this->Bp);

    float J = this->J;
    std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)> bondfunc = 
    [J](const Eigen::Vector2d &S1, const Eigen::Vector2d &S2) {
        return -J*S1.dot(S2);
    };

    if (!cluster_update)
        this->mut_mode = 0;

    Eigen::Vector3d v1; v1 << 1.,0.,0.;
    Eigen::Vector3d v2; v2 << 0.,1.,0.;
    this->add_bond(1,0,0,0,   v1, bondfunc);
    this->add_bond(-1,0,0,0, -v1, bondfunc);
    this->add_bond(0,1,0,0,   v2, bondfunc);
    this->add_bond(0,-1,0,0, -v2, bondfunc);
}

std::vector<double> SquareXYModel::vorticity() const {
    float v1 = 0;
    float v2 = 0;

    std::vector<std::vector<std::vector<float>>> phi = std::vector<std::vector<std::vector<float>>>(N,
                                                std::vector<std::vector<float>>(N,
                                                        std::vector<float>(L)));

    int i;
    for (int n1 = 0; n1 < N; n1++) {
        for (int n2 = 0; n2 < N; n2++) {
            for (int n3 = 0; n3 < L; n3++) {
                i = flat_idx(n1, n2, n3, 0);
                phi[n1][n2][n3] = 0.;
                phi[n1][n2][n3] = std::atan2(get_spin(i)[1], get_spin(i)[0]);
            }
        }
    }

    float p1; float p2; float p3; float p4;
    float w;
    for (int n1 = 0; n1 < N; n1++) {
        for (int n2 = 0; n2 < N; n2++) {
            for (int n3 = 0; n3 < L; n3++) {
                p1 = phi[n1][n2][n3]; p2 = phi[(n1+1)%N][n2][n3];
                p3 = phi[(n1+1)%N][(n2+1)%N][n3]; p4 = phi[n1][(n2+1)%N][n3];
                w = arg(exp(std::complex<float>(0., p2 - p1))) + arg(exp(std::complex<float>(0., p3 - p2)))
                    + arg(exp(std::complex<float>(0., p4 - p3))) + arg(exp(std::complex<float>(0., p1 - p4)));
                if (w > 0) { v1 += w; } else { v2 += w; }
            }
        }
    }
    
    return std::vector<double>{v1/(2*PI*N*N*L), v2/(2*PI*N*N*L)};
}

float SquareXYModel::p(int i) const {
    return std::atan2(get_spin(i)[1], get_spin(i)[0]);
}

float SquareXYModel::e1() const {
    float s = 0;
    for (int i = 0; i < V; i++) {
        s += std::cos(p(i) - p(neighbors[i][0].first));
    }
    return s/V;
}

float SquareXYModel::e2() const {
    float s = 0;
    for (int i = 0; i < V; i++) {
        s += std::sin(p(i) - p(neighbors[i][0].first));
    }
    return s/V;
}

float SquareXYModel::U2() const {
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


    Eigen::Vector2d S1;
    Eigen::Vector2d S2;
    int j;
    for (int i = 0; i < V; i++) {
        j = neighbors[i][0].first;

        S1 = get_spin(i);
        S2 = get_spin(j);

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
    float E = 0;

    // Onsite interactions
    E -= Bx*S[0] + By*S[1];
    return E;
}

void SquareXYModel::over_relaxation_mutation() {
    Eigen::Vector2d H; H << 0., 0.;
    int j;
    for (int n = 0; n < bonds.size(); n++) {
        j = neighbors[mut.i][n].first;
        H -= J*get_spin(j);
    }

    this->mut.dS = -2*get_spin(mut.i) + 2.*get_spin(mut.i).dot(H)/std::pow(H.norm(),2) * H;
}

void SquareXYModel::generate_mutation() {
    if (cluster_update)
        cluster_mutation(); 
    else {
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

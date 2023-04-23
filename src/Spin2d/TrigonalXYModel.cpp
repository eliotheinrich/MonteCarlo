#include "TrigonalXYModel.h"
#include <iostream>
#include <functional>
#include <string>

TrigonalXYModel::TrigonalXYModel(Params &params) : Spin2DModel(params) {
    this->N = params.get<int>("system_size");
    this->L = params.get<int>("layers", DEFAULT_LAYERS);
    this->J = params.get<float>("J");
    this->A = params.get<float>("A");

    Spin2DModel::init_params(1, N, N, L);

    float J = this->J;
    std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)> dotfunc = 
    [J](const Eigen::Vector2d &S1, const Eigen::Vector2d &S2) {
        return -J*S1.dot(S2);
    };

    Eigen::Vector3d v1; v1 << 1., 0., 0.;
    Eigen::Vector3d v2; v2 << 0.5, std::sqrt(3)/2., 0.;
    Eigen::Vector3d v3; v3 << 0.5, -std::sqrt(3)/2., 0.;
    this->add_bond(1, 0,0,0,v1, dotfunc);
    this->add_bond(-1,0,0,0,-v1,dotfunc);
    this->add_bond(0, 1,0,0,v2, dotfunc);
    this->add_bond(0,-1,0,0,-v2,dotfunc);
    this->add_bond(1,-1,0,0,v3, dotfunc);
    this->add_bond(-1,1,0,0,-v3,dotfunc);
    this->mut_mode = 0;
}

std::vector<double> TrigonalXYModel::vorticity() const {
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
                phi[n1][n2][n3] = atan2(get_spin(i)[1], get_spin(i)[0]);
            }
        }
    }

    float p1; float p2; float p3;
    float w;
    for (int n1 = 0; n1 < N; n1++) {
        for (int n2 = 0; n2 < N; n2++) {
            for (int n3 = 0; n3 < L; n3++) {
                p1 = phi[n1][n2][n3]; p2 = phi[(n1+1)%N][n2][n3];
                p3 = phi[(n1+1)%N][(n2+1)%N][n3];
                w = std::arg(exp(std::complex<float>(0., p2 - p1))) + std::arg(exp(std::complex<float>(0., p3 - p2)))
                    + std::arg(exp(std::complex<float>(0., p1 - p3)));
                if (w > 0) { v1 += w; } else { v2 += w; }
            }
        }
    }
    
    return std::vector<double>{v1/(2*PI*N*N*L), v2/(2*PI*N*N*L)};
}

double TrigonalXYModel::onsite_func(const Eigen::Vector2d &S) const {
    float phi = atan2(S[1], S[0]);
    return A*std::cos(6*phi);
}

void TrigonalXYModel::over_relaxation_mutation() {
    Eigen::Vector2d H; H << 0., 0.;
    int j;
    for (int n = 0; n < bonds.size(); n++) {
        j = neighbors[mut.i][n].first;
        H -= J*get_spin(j);
    }

    this->mut.dS = -2*get_spin(mut.i) + 2.*get_spin(mut.i).dot(H)/std::pow(H.norm(),2) * H;
}

void TrigonalXYModel::generate_mutation() {
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

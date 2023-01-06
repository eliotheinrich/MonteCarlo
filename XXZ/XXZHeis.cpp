#ifndef XXZ_
#define XXZ_

#include "XXZHeis.h"


XXZHeis::XXZHeis(int N, int L, float J, float K) : SpinModel(1, N, N, L) {
    this->N = N;
    this->L = L;
    this->J = J;
    this->K = K;


    std::function<float(Eigen::Vector3d, Eigen::Vector3d)> bondfunc = [J, K](Eigen::Vector3d S1, Eigen::Vector3d S2) {
        return -J*S1.dot(S2) + K*S1[2]*S2[2];
    };



    Eigen::Vector3d v1; v1 << 1.,0.,0.;
    Eigen::Vector3d v2; v2 << 0.,1.,0.;
    this->add_bond(1,0,0,0,   v1, bondfunc);
    this->add_bond(-1,0,0,0, -v1, bondfunc);
    this->add_bond(0,1,0,0,   v2, bondfunc);
    this->add_bond(0,-1,0,0, -v2, bondfunc);
}

XXZHeis* XXZHeis::clone() {
    XXZHeis* new_model = new XXZHeis(N, L, J, K);
    for (int i = 0; i < V; i++) {
        new_model->spins[i] = this->spins[i];
    }
    return new_model;
}

inline std::vector<double> XXZHeis::vorticity() {
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
                phi[n1][n2][n3] = atan2(spins[i][1], spins[i][0]);
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

double XXZHeis::onsite_func(const Eigen::Vector3d &S) const {
    // Onsite interactions
    return 0.;
}

#endif
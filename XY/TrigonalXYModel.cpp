#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../XYModel.cpp"
#include "../Utility.cpp"


using namespace std;
using namespace Eigen;

class TrigonalXYModel : public XYModel {
    public:
        int N;
        int L;
        float J;
        float A;

        int mut_mode;

        TrigonalXYModel(int N, int L, float J, float A) : XYModel(1, N, N, L) {
            this->N = N;
            this->L = L;
            this->J = J;
            this->A = A;

            std::function<float(Eigen::Vector2f, Eigen::Vector2f)> dotfunc = 
            [J](Eigen::Vector2f S1, Eigen::Vector2f S2) {
                return -J*S1.dot(S2);
            };

            Eigen::Vector3f v1; v1 << 1., 0., 0.;
            Eigen::Vector3f v2; v2 << 0.5, sqrt(3)/2., 0.;
            Eigen::Vector3f v3; v3 << 0.5, -sqrt(3)/2., 0.;
            this->add_bond(1, 0,0,0,v1, dotfunc);
            this->add_bond(-1,0,0,0,-v1,dotfunc);
            this->add_bond(0, 1,0,0,v2, dotfunc);
            this->add_bond(0,-1,0,0,-v2,dotfunc);
            this->add_bond(1,-1,0,0,v3, dotfunc);
            this->add_bond(-1,1,0,0,-v3,dotfunc);
            this->mut_mode = 0;
        }


        TrigonalXYModel* clone() {
            TrigonalXYModel* new_model = new TrigonalXYModel(N, L, J, A);
            for (int i = 0; i < V; i++) {
                new_model->spins[i] = this->spins[i];
            }
            return new_model;
        }

        inline std::vector<double> vorticity() {
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

            float p1; float p2; float p3;
            float w;
            for (int n1 = 0; n1 < N; n1++) {
                for (int n2 = 0; n2 < N; n2++) {
                    for (int n3 = 0; n3 < L; n3++) {
                        p1 = phi[n1][n2][n3]; p2 = phi[(n1+1)%N][n2][n3];
                        p3 = phi[(n1+1)%N][(n2+1)%N][n3];
                        w = arg(exp(complex<float>(0., p2 - p1))) + arg(exp(complex<float>(0., p3 - p2)))
                          + arg(exp(complex<float>(0., p1 - p3)));
                        if (w > 0) { v1 += w; } else { v2 += w; }
                    }
                }
            }
            
            return std::vector<double>{v1/(2*PI*N*N*L), v2/(2*PI*N*N*L)};
        }

        const float onsite_func(const Vector2f &S) {
            float phi = atan2(S[1], S[0]);
            return A*std::cos(6*phi);
        }

#ifndef CLUSTER_UPDATE
        void over_relaxation_mutation() {
            Eigen::Vector2f H; H << 0., 0.;
            int j;
            for (int n = 0; n < bonds.size(); n++) {
                j = neighbors[mut.i][n];
                H -= J*spins[j];
            }

            this->mut.dS = -2*spins[mut.i] + 2.*spins[mut.i].dot(H)/std::pow(H.norm(),2) * H;
        }
#endif

        void generate_mutation() {
#ifdef CLUSTER_UPDATE
            cluster_update(); 
#else
            mut.i = r() % V;

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
#endif
        }
};


#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include <complex>
#include "../XYModel.cpp"
#include "../Utility.cpp"


class SquareXYModel : public XYModel {
    public:
        int N;
        int L;
        int mut_mode;
        float J;
        float B;
        float Bp;
        float Bx;
        float By;

        SquareXYModel(int N, int L, float J, float B, float Bp) : XYModel(1, N, N, L) {
            this->N = N;
            this->L = L;
            this->J = J;
            this->B = B;
            this->Bp = Bp;
            this->Bx = B*cos(Bp);
            this->By = B*sin(Bp);

            std::function<float(Eigen::Vector2f, Eigen::Vector2f)> bondfunc = 
            [J](Eigen::Vector2f S1, Eigen::Vector2f S2) {
                return -J*S1.dot(S2);
            };

            this->mut_mode = 0;

            Eigen::Vector3f v1; v1 << 1.,0.,0.;
            Eigen::Vector3f v2; v2 << 0.,1.,0.;
            this->add_bond(1,0,0,0,   v1, bondfunc);
            this->add_bond(-1,0,0,0, -v1, bondfunc);
            this->add_bond(0,1,0,0,   v2, bondfunc);
            this->add_bond(0,-1,0,0, -v2, bondfunc);
        }

        SquareXYModel* clone() {
            SquareXYModel* new_model = new SquareXYModel(N, L, J, B, Bp);
            new_model->cluster = this->cluster;
            for (int i = 0; i < V; i++) {
                new_model->spins[i]  = spins[i];
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
                        phi[n1][n2][n3] = 0.;
                        phi[n1][n2][n3] = std::atan2(spins[i][1], spins[i][0]);
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
                        w = arg(exp(complex<float>(0., p2 - p1))) + arg(exp(complex<float>(0., p3 - p2)))
                          + arg(exp(complex<float>(0., p4 - p3))) + arg(exp(complex<float>(0., p1 - p4)));
                        if (w > 0) { v1 += w; } else { v2 += w; }
                    }
                }
            }
            
            return std::vector<double>{v1/(2*PI*N*N*L), v2/(2*PI*N*N*L)};
        }

        float p(int i) {
            return std::atan2(spins[i][1], spins[i][0]);
        }

        float e1() {
            float s = 0;
            for (int i = 0; i < V; i++) {
                s += std::cos(p(i) - p(neighbors[i][0]));
            }
            return s/V;
        }

        float e2() {
            float s = 0;
            for (int i = 0; i < V; i++) {
                s += std::sin(p(i) - p(neighbors[i][0]));
            }
            return s/V;
        }

        float U2() {
            return e1() - V/T*pow(e2(), 2);
        }

        inline std::vector<double> twist_stiffness() {
            auto Y = XYModel::twist_stiffness();
            return std::vector<double>{Y[0], Y[1], Y[2], Y[3], U2(), e1(), pow(e2(),4)};
        }

        const float onsite_func(const Eigen::Vector2f &S) {
            float E = 0;

            // Onsite interactions
            E -= Bx*S[0] + By*S[1];
            return E;
        }

        void over_relaxation_mutation() {
            Eigen::Vector2f H; H << 0., 0.;
            int j;
            for (int n = 0; n < bonds.size(); n++) {
                j = neighbors[mut.i][n];
                H -= J*spins[j];
            }

            this->mut.dS = -2*spins[mut.i] + 2.*spins[mut.i].dot(H)/std::pow(H.norm(),2) * H;
        }

        void generate_mutation() {
            if (cluster) { 
                cluster_update(); 
            } else {
                mut.i = (mut.i + 1) % V;

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
};


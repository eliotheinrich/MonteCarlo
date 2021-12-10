#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include <complex>
#include "../XYModel.cpp"
#include "../Utility.cpp"


using namespace std;
using namespace Eigen;

class SquareXYModel : public XYModel {
    public:
        int N;
        int L;
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

            function<float(Vector2f, Vector2f)> bondfunc_J = [J](Vector2f S1, Vector2f S2) {
                return -J*S1.dot(S2);
            };



            Vector3f v1; v1 << 1.,0.,0.;
            Vector3f v2; v2 << 0.,1.,0.;
            this->add_bond(Bond{1,0,0,0,v1,bondfunc_J});
            this->add_bond(Bond{-1,0,0,0,-v1,bondfunc_J});
            this->add_bond(Bond{0,1,0,0,v2,bondfunc_J});
            this->add_bond(Bond{0,-1,0,0,-v2,bondfunc_J});
        }

        SquareXYModel* clone() {
            SquareXYModel* new_model = new SquareXYModel(N, L, J, B, Bp);
            new_model->random_selection = random_selection;
            for (int n1 = 0; n1 < N; n1++) {
                for (int n2 = 0; n2 < N; n2++) {
                    for (int n3 = 0; n3 < L; n3++) {
                        new_model->spins[n1][n2][n3][0] = this->spins[n1][n2][n3][0];
                    }
                }
            }
            return new_model;
        }

        inline vector<double> vorticity() {
            float v1 = 0;
            float v2 = 0;

            vector<vector<vector<float>>> phi = vector<vector<vector<float>>>(N,
                                                        vector<vector<float>>(N,
                                                                vector<float>(L)));

            for (int n1 = 0; n1 < N; n1++) {
                for (int n2 = 0; n2 < N; n2++) {
                    for (int n3 = 0; n3 < L; n3++) {
                        phi[n1][n2][n3] = 0.;
                        phi[n1][n2][n3] = atan2(spins[n1][n2][n3][0][1], spins[n1][n2][n3][0][0]);
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
            
            return vector<double>{v1/(2*PI*N*N*L), v2/(2*PI*N*N*L)};
        }

        const float onsite_energy(int n1, int n2, int n3, int s) {
            float E = 0;

            // Onsite interactions
            //E -= this->Bx*this->spins[n1][n2][n3][s][0] + this->By*this->spins[n1][n2][n3][s][1];

            return E;
        }

        const float bond_energy(int n1, int n2, int n3, int s) {
            float E = 0;

            // NN interactions

            // PBC boundary conditions
            E -= this->spins[n1][n2][n3][s].dot(this->spins[(n1+1)%N][n2][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[mod(n1-1,N)][n2][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[n1][(n2+1)%N][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[n1][mod(n2-1,N)][n3][s]);

            return 0.5*J*E;
        }
};


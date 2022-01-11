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

            function<float(Vector2f, Vector2f)> bondfunc = [J](Vector2f S1, Vector2f S2) {
                return -J*S1.dot(S2);
            };



            Vector3f v1; v1 << 1.,0.,0.;
            Vector3f v2; v2 << 0.,1.,0.;
            this->add_bond(Bond{1,0,0,0,   v1, bondfunc});
            this->add_bond(Bond{-1,0,0,0, -v1, bondfunc});
            this->add_bond(Bond{0,1,0,0,   v2, bondfunc});
            this->add_bond(Bond{0,-1,0,0, -v2, bondfunc});
        }

        SquareXYModel* clone() {
            SquareXYModel* new_model = new SquareXYModel(N, L, J, B, Bp);
            for (int i = 0; i < V; i++) {
                new_model->spins[i]  = spins[i];
            }
            return new_model;
        }

        inline vector<double> vorticity() {
            float v1 = 0;
            float v2 = 0;

            vector<vector<vector<float>>> phi = vector<vector<vector<float>>>(N,
                                                        vector<vector<float>>(N,
                                                                vector<float>(L)));

            int i;
            for (int n1 = 0; n1 < N; n1++) {
                for (int n2 = 0; n2 < N; n2++) {
                    for (int n3 = 0; n3 < L; n3++) {
                        i = flat_idx(n1, n2, n3, 0);
                        phi[n1][n2][n3] = 0.;
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
                        w = arg(exp(complex<float>(0., p2 - p1))) + arg(exp(complex<float>(0., p3 - p2)))
                          + arg(exp(complex<float>(0., p4 - p3))) + arg(exp(complex<float>(0., p1 - p4)));
                        if (w > 0) { v1 += w; } else { v2 += w; }
                    }
                }
            }
            
            return vector<double>{v1/(2*PI*N*N*L), v2/(2*PI*N*N*L)};
        }

        const float onsite_energy(int i) {
            float E = 0;

            // Onsite interactions
            E -= Bx*spins[i][0] + By*spins[i][1];
            return E;
        }
};


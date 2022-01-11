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

        TrigonalXYModel(int N, int L, float J, float A) : XYModel(1, N, N, L) {
            this->N = N;
            this->L = L;
            this->J = J;
            this->A = A;

            function<float(Vector2f, Vector2f)> dotfunc = [J](Vector2f S1, Vector2f S2) {
                return -J*S1.dot(S2);
            };

            Vector3f v1; v1 << 1., 0., 0.;
            Vector3f v2; v2 << 0.5, sqrt(3)/2., 0.;
            Vector3f v3; v3 << 0.5, -sqrt(3)/2., 0.;
            this->add_bond(Bond{1,0,0,0,v1,dotfunc});
            this->add_bond(Bond{-1,0,0,0,-v1,dotfunc});
            this->add_bond(Bond{0,1,0,0,v2,dotfunc});
            this->add_bond(Bond{0,-1,0,0,-v2,dotfunc});
            this->add_bond(Bond{1,-1,0,0,v3,dotfunc});
            this->add_bond(Bond{-1,1,0,0,-v3,dotfunc});
        }


        TrigonalXYModel* clone() {
            TrigonalXYModel* new_model = new TrigonalXYModel(N, L, J, A);
            for (int i = 0; i < V; i++) {
                new_model->spins[i] = this->spins[i];
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
            
            return vector<double>{v1/(2*PI*N*N*L), v2/(2*PI*N*N*L)};
        }

        const float onsite_energy(int i) {
            float phi = atan2(spins[i][1], spins[i][0]);
            return A*cos(6*phi);
        }
};


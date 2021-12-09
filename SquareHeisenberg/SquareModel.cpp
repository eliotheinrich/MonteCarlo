#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../SpinModel.cpp"
#include "../Utility.cpp"


using namespace std;
using namespace Eigen;


class SquareModel : public SpinModel {
    public:
        int N;
        int L;
        float J;
        float A;


        SquareModel(int N, int L, float J, float A) : SpinModel(1, N, N, L) {

            this->N = N;
            this->L = L;
            this->J = J;
            this->A = A;
            function<float(Vector3f S1, Vector3f S2)> bondfunc = [J](Vector3f S1, Vector3f S2) {
                return -J*S1.dot(S2);
            };

            Vector3f v1; v1 << 1., 0., 0.;
            Vector3f v2; v2 << 0., 1., 0.;
            this->add_bond(Bond{1,0,0,0, v1, bondfunc});
            this->add_bond(Bond{-1,0,0,0, -v1, bondfunc});
            this->add_bond(Bond{0,1,0,0, v2, bondfunc});
            this->add_bond(Bond{0,-1,0,0, -v2, bondfunc});
        }

        SquareModel* clone() {
            SquareModel* new_model = new SquareModel(N, L, J, A);
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

        const float onsite_energy(int n1, int n2, int n3, int s) {
            return A*pow(this->spins[n1][n2][n3][s][2], 2);
        }

        const float bond_energy(int n1, int n2, int n3, int s) {
            float E = 0;

            // NN interactions
            E -= spins[n1][n2][n3][s].dot(spins[n1][(n2+1)%N][n3][s]);
            E -= spins[n1][n2][n3][s].dot(spins[n1][mod(n2-1, N)][n3][s]);
            E -= spins[n1][n2][n3][s].dot(spins[(n1+1)%N][n2][n3][s]);
            E -= spins[n1][n2][n3][s].dot(spins[mod(n1-1, N)][n2][n3][s]);

            return 0.5*J*E;
        }
};


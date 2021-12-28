#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../IsingModel.cpp"
#include "../Utility.cpp"


using namespace std;
using namespace Eigen;


class SquareIsingModel : public IsingModel {
    public:
        int N;
        int L;
        float J;
        float B;

        SquareIsingModel(int N, int L, float J, float B) : IsingModel(N, N, L) {
            this->N = N;
            this->L = L;
            this->J = J;
            this->B = B;
        }

        SquareIsingModel *clone() {
            SquareIsingModel *new_model = new SquareIsingModel(N, L, J, B);
            for (int n1 = 0; n1 < N; n1++) {
                for (int n2 = 0; n2 < N; n2++) {
                    for (int n3 = 0; n3 < L; n3++) {
                        new_model->spins[n1][n2][n3] = this->spins[n1][n2][n3];
                    }
                }
            }
            return new_model;
        }

        const float onsite_energy(int n1, int n2, int n3) {
            float E = 0;

            // Onsite interactions
            E -= B*this->spins[n1][n2][n3];

            return E;
        }

        const float bond_energy(int n1, int n2, int n3) {
            float E = 0;

            // NN interactions

            // PBC boundary conditions
            E -= J*this->spins[n1][n2][n3]*this->spins[(n1+1)%N][n2][n3];
            E -= J*this->spins[n1][n2][n3]*this->spins[mod(n1-1,N)][n2][n3];

            E -= J*this->spins[n1][n2][n3]*this->spins[n1][(n2+1)%N][n3];
            E -= J*this->spins[n1][n2][n3]*this->spins[n1][mod(n2-1,N)][n3];


            return 0.5*E;
        }

};


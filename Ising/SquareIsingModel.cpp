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

//            E -= J*this->spins[n1][n2][n3]*this->spins[n1][n2][(n3+1)%L];
//            E -= J*this->spins[n1][n2][n3]*this->spins[n1][n2][mod(n3-1,L)];

            return 0.5*E;
        }

};


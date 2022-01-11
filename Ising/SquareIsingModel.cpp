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
            for (int i = 0; i < V; i++) {
                new_model->spins[i] = this->spins[i];
            }
            return new_model;
        }

        const float onsite_energy(int i) {
            float E = 0;

            // Onsite interactions
            E -= B*this->spins[i];

            return E;
        }

        const float bond_energy(int i) {
            float E = 0;

            Vector3i idxs = tensor_idx(i);
            int n1 = idxs[0]; int n2 = idxs[1]; int n3 = idxs[2]; 

            // NN interactions
            E -= J*spins[i]*spins[flat_idx(mod(n1+1, N), n2, n3)];
            E -= J*spins[i]*spins[flat_idx(mod(n1-1, N), n2, n3)];

            E -= J*spins[i]*spins[flat_idx(n1, mod(n2+1, N), n3)];
            E -= J*spins[i]*spins[flat_idx(n1, mod(n2-1, N), n3)];

            return 0.5*E;
        }

};


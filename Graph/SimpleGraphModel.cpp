#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../GraphModel.cpp"
#include "../Utility.cpp"

class SimpleGraphModel : public GraphModel {
    public:
        int N;

        SimpleGraphModel(int N, float J) : GraphModel(N) {
            this->N = N;
            this->J = J;
        }

        SimpleGraphModel *clone() {
            SquareIsingModel *new_model = new SimpleGraphModel(N, J);
        }

        const float onsite_energy(int i) {
            float E = 0;

            return E;
        }

        const float bond_energy(int i) {
            float E = 0;
            E += deg(i);
            return J*E;
        }

};


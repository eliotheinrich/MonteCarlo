#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include "../GraphModel.cpp"
#include "../Utility.cpp"

class SimpleGraphModel : public GraphModel {
    public:
        int N;
        float J;

        SimpleGraphModel(int N, float J) : GraphModel(N) {
            this->N = N;
            this->J = J;
        }

        SimpleGraphModel *clone() {
            SimpleGraphModel *new_model = new SimpleGraphModel(N, J);
            return new_model;
        }

        const double onsite_energy(int i) {
            float E = 0;

            return E;
        }

        const double bond_energy(int i) {
            float E = 0;
            E += deg(i);
            return -J*E;
        }

};


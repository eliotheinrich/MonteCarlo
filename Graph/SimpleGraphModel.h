#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include "../GraphModel.cpp"
#include "../Utility.cpp"

class SimpleGraphModel : public GraphModel {
    public:
        int N;

        SimpleGraphModel(int N, float J);

        SimpleGraphModel *clone();

        const float onsite_energy(int i);
        const float bond_energy(int i);
};


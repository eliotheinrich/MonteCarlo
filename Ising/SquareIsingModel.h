#ifndef SQUAREISING_H
#define SQUAREISING_H

#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../IsingModel.h"
#include "../Utility.cpp"

class SquareIsingModel : public IsingModel {
    public:
        int N;
        int L;
        float J;
        float B;

        SquareIsingModel(int N, int L, float J, float B);

        SquareIsingModel *clone();

        double onsite_energy(int i) const;

        double bond_energy(int i) const;
};

#endif
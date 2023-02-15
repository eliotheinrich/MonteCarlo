#ifndef SQUAREISING_H
#define SQUAREISING_H

#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "IsingModel.h"
#include "Utility.h"

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


        virtual std::map<std::string, int> get_int_params() const;
        virtual std::map<std::string, double> get_double_params() const;
};

#endif
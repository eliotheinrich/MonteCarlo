#ifndef SQUARECLOCK_H
#define SQUARECLOCK_H

#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include <complex>
#include "../ClockModel.h"
#include "../Utility.cpp"

template <int q>
class SquareClockModel : public ClockModel<q> {
    public:
        int N;
        int L;
        float J;

        float bond_table[q][q];

        SquareClockModel(int N, int L, float J);

        SquareClockModel<q>* clone();

        virtual double onsite_energy(int i) const;
};

#include "SquareClockModel.cpp"

#endif
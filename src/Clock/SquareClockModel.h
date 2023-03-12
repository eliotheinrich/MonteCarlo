#ifndef SQUARECLOCK_H
#define SQUARECLOCK_H

#include <Eigen/Dense>
#include "ClockModel.h"

#define DEFAULT_LAYERS 1

template <int q>
class SquareClockModel : public ClockModel<q> {
    public:
        int N;
        int L;
        float J;

        float bond_table[q][q];

        SquareClockModel(Params &params);

        virtual MCModel* clone(Params &params) { return new SquareClockModel<q>(params); };

        virtual double onsite_energy(int i) const;
};

#include "SquareClockModel.cpp"

#endif
#ifndef SQUAREISING_H
#define SQUAREISING_H

#include "IsingModel.h"

#define DEFAULT_LAYERS 1

class SquareIsingModel : public IsingModel {
    public:
        int N;
        int L;
        float J;
        float B;

        SquareIsingModel(Params &params);

        virtual double onsite_energy(int i) const;
        virtual double bond_energy(int i) const;

        CLONE(MCModel, SquareIsingModel)
};

#endif
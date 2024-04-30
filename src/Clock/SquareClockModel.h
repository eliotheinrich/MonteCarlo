#pragma once

#include <Eigen/Dense>
#include "ClockModel.h"

#define DEFAULT_LAYERS 1

template <uint32_t q>
class SquareClockModel : public ClockModel<q> {
  public:
    uint32_t N;
    uint32_t L;
    double J;

    double bond_table[q][q];

    SquareClockModel(dataframe::Params &params, uint32_t num_threads);

    virtual double onsite_energy(uint32_t i) const override;
};

#include "SquareClockModel.cpp"


#pragma once

#include "IsingModel.h"

#define DEFAULT_LAYERS 1

class SquareIsingModel : public IsingModel {
  public:
    SquareIsingModel(dataframe::Params &params, uint32_t num_threads);

    virtual double onsite_energy(uint32_t i) const override;
    virtual double bond_energy(uint32_t i) const override;

  private:
    uint32_t N;
    uint32_t L;
    double J;
    double B;
};


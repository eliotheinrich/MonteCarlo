#pragma once

#include "MultibodyIsingModel.h"

class LDPCIsingModel : public MultibodyIsingModel {
  public:
    LDPCIsingModel(dataframe::ExperimentParams &params, uint32_t num_threads);

    virtual double onsite_energy(uint32_t i) const override;

  private:
    size_t L;
    double J;
    double impurity;
};


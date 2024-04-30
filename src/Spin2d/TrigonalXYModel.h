#pragma once

#include <vector>
#include "Spin2DModel.h"

#define DEFAULT_LAYERS 1

class TrigonalXYModel : public Spin2DModel {
  public:

    TrigonalXYModel(dataframe::Params &params, uint32_t num_threads);

    std::vector<double> vorticity() const;

    virtual double onsite_func(const Eigen::Vector2d &S) const override;

    void over_relaxation_mutation();
    virtual void generate_mutation() override;

  private:
    uint32_t N;
    uint32_t L;
    double J;
    double A;

    int mut_mode;
};


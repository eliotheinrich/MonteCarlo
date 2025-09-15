#pragma once

#include <vector>
#include <Eigen/Dense>
#include "Spin3DModel.h"

#define DEFAULT_LAYERS 1

class XXZHeis : public Spin3DModel {
  public:
    XXZHeis(dataframe::ExperimentParams &params, uint32_t num_threads);

    std::vector<double> vorticity() const;

    virtual double onsite_func(const Eigen::Vector3d &S) const override;

    virtual dataframe::SampleMap take_samples() const override;

  private: 
    uint32_t N;
    uint32_t L;
    double J;
    double K;
    double A;

    bool sample_helicity;
    bool sample_structure_factor;

    void add_structure_factor_samples(dataframe::SampleMap& samples) const;
};

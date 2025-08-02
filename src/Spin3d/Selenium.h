#pragma once

#include <vector>
#include "Spin3DModel.h"

class SeleniumModel : public Spin3DModel {
  public:
    SeleniumModel(dataframe::ExperimentParams &params, uint32_t num_threads);

    void over_relaxation_mutation();
    virtual void generate_mutation() override;

    virtual double onsite_func(const Eigen::Vector3d &S) const override;

    void add_sublattice_magnetization_samples(dataframe::SampleMap& samples) const;

    void add_structure_factor_samples(dataframe::SampleMap& samples) const;
    void add_staggered_structure_factor_samples(dataframe::SampleMap& samples) const;

    virtual dataframe::SampleMap take_samples() override;

  private:
    uint32_t N;

    double J1;
    double J2;
    double J2p;
    double K;
    double D1;
    double D2;

    bool sample_sublattice_magnetization;
    bool sample_structure_factor;
    bool sample_staggered_structure_factor;

    uint32_t mut_counter;
    uint32_t mut_mode;
};

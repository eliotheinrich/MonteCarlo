#pragma once

#include <vector>
#include "Spin3DModel.h"

class TrigonalModel : public Spin3DModel {
  public:
    TrigonalModel(dataframe::ExperimentParams &params, uint32_t num_threads);

    void over_relaxation_mutation();
    virtual void generate_mutation() override;

    virtual double onsite_func(const Eigen::Vector3d &S) const override;

    // --- FOR DYNAMIC UPDATES --- //
    void rotate_spin(uint32_t i, Eigen::Vector3d v, double p);
    Eigen::Vector3d molecular_field(uint32_t i) const;
    void dynamic_step(double dt);
    // ------ //

    // For computing structure factor
    double intensity(Eigen::Vector3d Q) const;

    void add_intensityx_samples(dataframe::SampleMap &samples) const;
    void add_intensityy_samples(dataframe::SampleMap &samples) const;
    void add_intensityz_samples(dataframe::SampleMap &samples) const;
    void add_layer_magnetization_samples(dataframe::SampleMap &samples) const;

    virtual dataframe::SampleMap take_samples() const override;

  private:
    uint32_t N;
    uint32_t L;
    double J1;
    double J2;
    double J3;
    double K1;
    double K2;
    double K3;
    Eigen::Vector3d B; 

    uint32_t fm_layers;

    Eigen::Matrix3d R;
    uint32_t mut_counter;
    uint32_t mut_mode;

    bool sample_magnetization;

    bool sample_helicity;

    bool sample_layer_magnetization;
    bool sample_intensityx;
    bool sample_intensityy;
    bool sample_intensityz;
    double max_L;
    double min_L;
    uint32_t intensity_resolution;
};

#pragma once

#include <vector>
#include "Spin3DModel.h"

class AltermagnetModel : public Spin3DModel {
  public:
    AltermagnetModel(dataframe::ExperimentParams &params, uint32_t num_threads);
    virtual void init() override;

    void over_relaxation_mutation();

    virtual double onsite_func(const Eigen::Vector3d &S) const override;
    virtual void annealing_callback(int epoch, int num_epochs) override {
      if (anneal) {
        double alpha = static_cast<double>(epoch)/num_epochs;
        B   =   B_i + alpha * (B_i - B_f);
        K   =   K_i + alpha * (K_i - K_f);
        J1  =  J1_i + alpha * (J1_i - J1_f);
        J2  =  J2_i + alpha * (J2_i - J2_f);
        J2p = J2p_i + alpha * (J2p_i - J2p_f);
        D1  =  D1_i + alpha * (D1_i - D1_f);
        D2  =  D2_i + alpha * (D2_i - D2_f);
      }
    }

    void add_sublattice_magnetization_samples(dataframe::SampleMap& samples) const;
    void add_structure_factor_samples(dataframe::SampleMap& samples) const;

    virtual dataframe::SampleMap take_samples() const override;

  private:
    uint32_t N;

    double J1;
    double J2;
    double J2p;
    double K;
    double D1;
    double D2;
    double B;

    double J1_i;
    double J2_i;
    double J2p_i;
    double K_i;
    double D1_i;
    double D2_i;
    double B_i;

    double J1_f;
    double J2_f;
    double J2p_f;
    double K_f;
    double D1_f;
    double D2_f;
    double B_f;

    bool sample_sublattice_magnetization;
    bool sample_structure_factor;

    bool anneal;

    uint32_t mut_counter;
    uint32_t mut_mode;

    std::string initial_state;
};

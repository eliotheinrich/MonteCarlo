#pragma once

#include <vector>
#include "Spin3DModel.h"

class AltermagnetModel : public Spin3DModel {
  public:
    AltermagnetModel(Params &params, uint32_t num_threads);
    virtual void init() override;

    void over_relaxation_mutation();

    virtual double onsite_func(const Eigen::Vector3d &S) const override;
    virtual void annealing_callback(int epoch, int num_epochs) override {
      if (anneal) {
        double alpha = static_cast<double>(epoch)/num_epochs;
        B   =         B_i + alpha * (B_f - B_i);
        K   =         K_i + alpha * (K_f - K_i);
        J1  =        J1_i + alpha * (J1_f - J1_i);
        J2  =        J2_i + alpha * (J2_f - J2_i);
        J2p =       J2p_i + alpha * (J2p_f - J2p_i);
        D1  =        D1_i + alpha * (D1_f - D1_i);
        D2  =        D2_i + alpha * (D2_f - D2_i);
        temperature = T_i + alpha * (T_f - T_i);
      }
    }

  private:
    uint32_t N;

    double T_i;
    double T_f;

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

    bool anneal;

    uint32_t mut_counter;
    uint32_t mut_mode;

    std::string initial_state;
};

#pragma once

#include <vector>
#include "Spin3DModel.h"

class TrigonalModel : public Spin3DModel {
  public:
    TrigonalModel(Params &params, uint32_t num_threads);

    void over_relaxation_mutation();
    virtual void generate_mutation() override;

    virtual double onsite_func(const Eigen::Vector3d &S) const override;

    // --- FOR DYNAMIC UPDATES --- //
    void rotate_spin(uint32_t i, Eigen::Vector3d v, double p);
    Eigen::Vector3d molecular_field(uint32_t i) const;
    void dynamic_step(double dt);
    // ------ //

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

};

#pragma once

#include <vector>
#include "Spin3DModel.h"

class HelixModel : public Spin3DModel {
  public:
    HelixModel(Params &params, uint32_t num_threads);

    virtual double onsite_func(const Eigen::Vector3d &S) const override {
      return 0.0;
    }

    virtual void annealing_callback(int epoch, int num_epochs) override {
      if (anneal) {
        double alpha = static_cast<double>(epoch)/num_epochs;
        temperature = max_temp - alpha * (max_temp- min_temp);
      }
    }

  private:
    uint32_t N;

    double J1;
    double J2;
    double D1;
    double D2;

    bool anneal;
    double min_temp;
    double max_temp;
};

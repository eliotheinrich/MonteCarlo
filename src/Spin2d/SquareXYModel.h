#pragma once

#include <vector>
#include <Eigen/Dense>
#include "Spin2DModel.h"

#define DEFAULT_LAYERS 1

class SquareXYModel : public Spin2DModel {
  public:
    SquareXYModel(dataframe::ExperimentParams &params, uint32_t num_threads);

    std::vector<std::pair<size_t, bool>> get_vortices() const;
    std::vector<double> vorticity() const;

    double p(uint32_t i) const;
    double e1() const;
    double e2() const;
    double U2() const;
    virtual std::vector<double> twist_stiffness() const override;

    virtual double onsite_func(const Eigen::Vector2d &S) const override;

    void over_relaxation_mutation();
    virtual void generate_mutation() override;

    virtual Texture get_texture() const override;

  private:
    uint32_t N;
    uint32_t L;
    int mut_mode;
    double J;
    double B;
    double Bp;
    double Bx;
    double By;
};

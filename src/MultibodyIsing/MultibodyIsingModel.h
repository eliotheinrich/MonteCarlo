#pragma once

#include "MonteCarlo.h"

#include <Eigen/Dense>

#include <vector>

class MultibodyIsingModel : public MonteCarloSimulator {
  // Generic 3D Ising model
  // Must be supplied with number of sublattices
  public:
    size_t L;
    uint64_t V;

    struct MultibodyIsingTerm {
      double J;
      std::vector<size_t> inds;
    };

    double acceptance;
    std::vector<int> spins;
    std::vector<double> onsite_potential;
    std::vector<std::vector<MultibodyIsingTerm>> terms;

    MultibodyIsingModel(dataframe::Params &params, uint32_t num_threads) : MonteCarloSimulator(params, num_threads) {}
    MultibodyIsingModel()=default;
    virtual ~MultibodyIsingModel()=default;

    void init(size_t L);
    void add_term(const std::vector<size_t>& inds, double J);

    size_t idx(size_t x, size_t y) const; 
    std::pair<int, int> coordinates(size_t i) const;

    void randomize_spins();

    double get_magnetization() const;
    virtual double onsite_energy(uint32_t i) const;
    virtual double bond_energy(uint32_t i) const;
    
    virtual uint64_t system_size() const override;
    virtual void generate_mutation() override;
    virtual void accept_mutation() override;
    virtual void reject_mutation() override;
    virtual double energy() const override;
    virtual double energy_change() override;

    virtual dataframe::data_t take_samples() override;

  protected:
    struct MultibodyIsingMutation {
      uint32_t i;
    };

    // Mutation being considered is stored as an attribute of the model
    MultibodyIsingMutation mut;
};


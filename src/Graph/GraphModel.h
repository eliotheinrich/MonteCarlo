#pragma once

#include "MonteCarlo.h"

#include <vector>

class GraphModel : public MonteCarloSimulator {
  // Generic graph model
  public:
    uint64_t N;

    double acceptance;
    std::vector<std::vector<uint32_t>> edges;
    std::vector<int> vals;

    GraphModel(dataframe::ExperimentParams &params, uint32_t num_threads);
    GraphModel()=default;
    virtual ~GraphModel()=default;

    void init(uint64_t N);
    virtual uint64_t system_size() const override {
      return N;
    }

    virtual void generate_mutation() override;
    virtual void accept_mutation() override;
    virtual void reject_mutation() override;

    virtual double onsite_energy(uint32_t i) const=0;
    virtual double bond_energy(uint32_t i) const=0;

    virtual double energy() const override;
    virtual double energy_change() override;

    void toggle_edge(uint32_t i, uint32_t j);
    uint32_t deg(uint32_t i) const;

    double get_connectivity() const;
    virtual dataframe::SampleMap take_samples() const override;

  protected:
    struct GraphMutation {
      uint32_t i;
      uint32_t j;
    };

    // Mutation being considered is stored as an attribute of the model
    GraphMutation mut;
};


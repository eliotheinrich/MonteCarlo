#pragma once

#include <MonteCarlo.h>

#include <Eigen/Dense>

#include <vector>

class IsingModel : public MonteCarloSimulator {
  // Generic 3D Ising model
  // Must be supplied with number of sublattices
  public:
    uint32_t N1;
    uint32_t N2;
    uint32_t N3;
    uint64_t V;

    double acceptance;
    std::vector<double> spins;

    IsingModel(dataframe::ExperimentParams &params, uint32_t num_threads) : MonteCarloSimulator(params, num_threads) {}
    IsingModel()=default;
    virtual ~IsingModel()=default;

    void init(uint32_t N1, uint32_t N2, uint32_t N3);

    uint32_t flat_idx(uint32_t n1, uint32_t n2, uint32_t n3) const;

    Eigen::Vector3i tensor_idx(uint32_t i) const;

    void randomize_spins();

    double get_magnetization() const;
    virtual double onsite_energy(uint32_t i) const { return 0.0; }
    virtual double bond_energy(uint32_t i) const { return 0.0; }
    
    virtual uint64_t system_size() const override;
    virtual void generate_mutation() override;
    virtual void accept_mutation() override;
    virtual void reject_mutation() override;
    virtual double energy() const override;
    virtual double energy_change() override;

    virtual dataframe::SampleMap take_samples() const override;

    void save_spins(const std::string& filename);

  protected:
    struct IsingMutation {
      uint32_t i;
    };

    // Mutation being considered is stored as an attribute of the model
    IsingMutation mut;
};


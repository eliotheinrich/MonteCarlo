#pragma once

#include <MonteCarlo.h>

#include "Lattice.hpp"

#include <unordered_set>
#include <stack>

using Spin2D = Eigen::Vector2d;

class Spin2DModel : public MonteCarloSimulator {
  // Generic 2d spin model in up to 3d lattice
  public:
    uint32_t V;

    Spin2DModel(dataframe::ExperimentParams &params, uint32_t num_threads);
    virtual ~Spin2DModel()=default;

    void init(const Lattice<Spin2D>& lattice);

    virtual uint64_t system_size() const override {
      if (cluster_update) {
        return 1;
      }

      return V;
    }

    void randomize_spins();

    Spin2D get_spin(uint32_t i) const {
      Spin2D s = lattice.spins[i];
      if (cluster_update) {
        s = s0.transpose() * s;
      }
      return s;
    }

    void set_spin(uint32_t i, const Spin2D& spin) {
      lattice.spins[i] = spin;
    }

    virtual std::vector<double> twist_stiffness() const;

    Eigen::Vector2d get_magnetization() const;

    virtual double onsite_energy(int i) const;
    virtual double bond_energy(int i) const;

    void metropolis_mutation();
    void cluster_mutation();

    virtual void generate_mutation() override;
    virtual void accept_mutation() override;
    virtual void reject_mutation() override;

    virtual double energy() const override;
    virtual double energy_change() override;

    virtual double onsite_func(const Eigen::Vector2d& S) const = 0;

    void add_magnetization_samples(dataframe::SampleMap &samples) const;
    void add_helicity_samples(dataframe::SampleMap &samples) const;
    void add_intensityx_samples(dataframe::SampleMap& samples) const;
    void add_intensityy_samples(dataframe::SampleMap& samples) const;
    void add_intensityz_samples(dataframe::SampleMap& samples) const;

    virtual dataframe::SampleMap take_samples() const override;

  protected:
    Lattice<Spin2D> lattice;

    // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
    // dS must conserve the norm of S[n1,n2,n3,s]
    struct Spin2DMutation {
      uint32_t i;
      Eigen::Vector2d dS;
    };

    static constexpr double alpha = 0.01;

    bool cluster_update;

    // Mutation being considered is stored as an attribute of the model
    Spin2DMutation mut;

  private:
    bool sample_energy;
    bool sample_magnetization;
    bool sample_helicity;

    bool sample_intensityx;
    bool sample_intensityy;
    bool sample_intensityz;
    double max_L;
    double min_L;
    uint32_t intensity_resolution;

    double acceptance;
    double sigma;

    std::unordered_set<uint32_t> s;
    Eigen::Matrix2d s0;
};

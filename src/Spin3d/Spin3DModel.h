#pragma once

#include "MonteCarlo.hpp"

#include "Lattice.hpp"

#include <vector>
#include <unordered_set>
#include <math.h>
#include <random>

using Spin3D = Eigen::Vector3d;

class GaussianDist {
  public:
    GaussianDist() {}
    GaussianDist(double mean, double std);

    double sample();

  private:
    std::minstd_rand rd;
    std::default_random_engine gen;
    std::normal_distribution<> dist;
};

class Spin3DModel : public MonteCarloSimulator {
  // Generic 3D Heisenberg model
  // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
  // dS must conserve the norm of S[n1,n2,n3,s]
  public:
    uint64_t V;

    static constexpr int METROPOLIS = 0;
    static constexpr int ADAPTIVE_METROPOLIS = 1;
    static constexpr int CLUSTER = 2;

    Spin3DModel(Params &params, uint32_t num_threads);
    virtual ~Spin3DModel()=default;

    void init(const Lattice<Spin3D>& lattice);
    virtual void init() { }

    // TODO change this
    virtual uint64_t system_size() const override {
      if (mutation_type == CLUSTER) {
        return 1;
      }

      return V;
    }

    void set_spin(uint32_t i, const Spin3D& S) { 
      lattice.spins[i] = S; 
    }

    Spin3D get_spin(uint32_t i) const { 
      Spin3D s = lattice.spins[i];
      if (mutation_type == CLUSTER) {
        s = s0.transpose() * s;
      }
      return s;
    }

    void randomize_spins();

    static std::vector<double> twist_terms(std::vector<double> dE);
    std::vector<double> twist_derivatives(uint32_t i) const;
    std::vector<double> twist_derivatives() const;

    Eigen::Vector3d get_magnetization() const;

    std::vector<double> correlation_function(uint32_t i, uint32_t a, uint32_t b) const;
    std::vector<double> full_correlation_function(uint32_t i) const;
    double skyrmion_density(uint32_t i) const;

    std::vector<double> skyrmion_correlation_function(uint32_t i) const;

    static Eigen::Vector3d mutate_spin3d(const Spin3D& spin, GaussianDist& dist, double sigma) {
      Eigen::Vector3d Gamma;
      Gamma << dist.sample(), dist.sample(), dist.sample();
      return (spin + sigma*Gamma).normalized();
    }

    void metropolis_mutation();
    void adaptive_metropolis_mutation();
    void cluster_mutation();

    virtual double energy() const override;
    virtual double energy_change() override;
    virtual void generate_mutation() override;
    virtual void accept_mutation() override;
    virtual void reject_mutation() override;

    virtual double onsite_func(const Spin3D& S) const=0;
    virtual double onsite_energy(uint32_t i) const;
    virtual double bond_energy(uint32_t i) const;

    LatticeGraph to_graph() const { 
      return lattice.to_graph();
    }

  protected:
    Lattice<Spin3D> lattice;

    Eigen::Matrix3d s0;

    struct Spin3DMutation {
      uint32_t i;
      Spin3D dS;
    };

    int mutation_type;

    // Mutation being considered is stored as an attribute of the model
    Spin3DMutation mut;

    // Internal normally distributed random number generator
    GaussianDist dist;

    static constexpr double alpha = 0.01;

  private:
    uint64_t nsteps;
    uint64_t accepted;
    double acceptance;
    double sigma;

    std::unordered_set<int> s;
};

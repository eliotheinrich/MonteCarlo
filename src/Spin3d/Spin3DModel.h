#pragma once

#include "MonteCarlo.h"

#include <Eigen/Dense>

#include <vector>
#include <unordered_set>
#include <math.h>
#include <random>

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

    Spin3DModel(dataframe::Params &params, uint32_t num_threads);
    virtual ~Spin3DModel()=default;

    void init(uint32_t sl, uint32_t N1, uint32_t N2, uint32_t N3);

    virtual uint64_t system_size() const override {
      if (cluster_update) {
        return 1;
      }

      return V;
    }

    void set_spin(uint32_t i, Eigen::Vector3d S) { 
      spins[i] = S; 
    }

    Eigen::Vector3d get_spin(uint32_t i) const { 
      return cluster_update ? s0.transpose()*spins[i] : spins[i]; 
    }

    void randomize_spins();

    uint32_t flat_idx(uint32_t n1, uint32_t n2, uint32_t n3, uint32_t s) const {
      return n1 + N1*(n2 + N2*(n3 + N3*s));
    }

    Eigen::Vector4i tensor_idx(uint32_t i) const {
      uint32_t n1 = i % N1;
      i = i / N1;
      uint32_t n2 = i % N2;
      i = i / N2;
      uint32_t n3 = i % N3;
      i = i / N3;
      uint32_t s = i % sl;
      Eigen::Vector4i v; v << n1, n2, n3, s;
      return v;
    }

    void add_bond(
      int d1, 
      int d2, 
      int d3, 
      int ds, 
      Eigen::Vector3d v, 
      std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc
    ) {
      add_bond(d1, d2, d3, ds, v, bondfunc, [](uint32_t, uint32_t){ return true; });
    }

    void add_bond(
      int d1, 
      int d2, 
      int d3, 
      int ds, 
      Eigen::Vector3d v, 
      std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc,
      std::function<bool(uint32_t, uint32_t)> bond_filter
    );

    static std::vector<double> twist_terms(std::vector<double> dE);
    std::vector<double> twist_derivatives(uint32_t i) const;
    std::vector<double> twist_derivatives() const;

    Eigen::Vector3d get_magnetization() const;

    std::vector<double> correlation_function(uint32_t i, uint32_t a, uint32_t b) const;
    std::vector<double> full_correlation_function(uint32_t i) const;
    double skyrmion_density(uint32_t i) const;

    std::vector<double> skyrmion_correlation_function(uint32_t i) const;

    void cluster_mutation();
    void metropolis_mutation();

    virtual double energy() const override;
    virtual double energy_change() override;
    virtual void generate_mutation() override;
    virtual void accept_mutation() override;
    virtual void reject_mutation() override;

    virtual double onsite_func(const Eigen::Vector3d& S) const=0;
    virtual double onsite_energy(uint32_t i) const;
    virtual double bond_energy(uint32_t i) const;


    // Saves current spin configuration
    void save_spins(const std::string& filename);
    bool load_spins(const std::string& filename);

    void add_magnetization_samples(dataframe::data_t &samples) const;
    void add_helicity_samples(dataframe::data_t &samples) const;
    virtual dataframe::data_t take_samples() override;

  protected:
    struct Spin3DMutation {
      uint32_t i;
      Eigen::Vector3d dS;
    };

    struct Spin3DBond {
      int d1;
      int d2;
      int d3;
      int ds;
      Eigen::Vector3d v;
      std::function<double(const Eigen::Vector3d&, const Eigen::Vector3d&)> bondfunc;
    };   

    bool cluster_update;

    // Mutation being considered is stored as an attribute of the model
    Spin3DMutation mut;

    std::vector<std::vector<Bond>> neighbors;
    std::vector<Spin3DBond> bonds;

    static constexpr double alpha = 0.01;
    std::vector<Eigen::Matrix3d> R1s;
    std::vector<Eigen::Matrix3d> R2s;
    std::vector<Eigen::Matrix3d> R3s;

  private:
    bool sample_energy;
    bool sample_magnetization;
    bool sample_helicity;

    uint32_t sl;
    uint32_t N1;
    uint32_t N2;
    uint32_t N3;

    BoundaryCondition bcx;
    BoundaryCondition bcy;
    BoundaryCondition bcz;

    uint64_t nsteps;
    uint64_t accepted;
    double acceptance;
    double sigma;
    std::vector<Eigen::Vector3d> spins;

    // Need to store bond filters passed to add_bond so that neighbor matrix can 
    // be later filled out bit init
    std::vector<std::function<bool(uint32_t, uint32_t)>> bond_filters;

    std::unordered_set<int> s;
    Eigen::Matrix3d s0;

    // Internal normally distributed random number generator
    GaussianDist dist;
};

#pragma once

#include <MonteCarlo.h>

#include "Lattice.hpp"

#include <vector>
#include <unordered_set>
#include <math.h>
#include <random>

#include <finufft.h>

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

std::pair<std::vector<double>, std::vector<double>> split_complex_output(const std::vector<std::complex<double>>& complex_data);
std::pair<std::vector<double>, std::vector<double>> fft2d_channel(std::vector<double>& x, std::vector<double>& y, std::vector<std::complex<double>>& input, int N, int s=1);

class Spin3DModel : public MonteCarloSimulator {
  // Generic 3D Heisenberg model
  // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
  // dS must conserve the norm of S[n1,n2,n3,s]
  public:
    uint64_t V;

    Spin3DModel(dataframe::ExperimentParams &params, uint32_t num_threads);
    virtual ~Spin3DModel()=default;

    void init(const Lattice<Spin3D>& lattice);

    // TODO change this
    virtual uint64_t system_size() const override {
      if (cluster_update) {
        return 1;
      }

      return V;
    }

    void set_spin(uint32_t i, const Spin3D& S) { 
      lattice.spins[i] = S; 
    }

    Spin3D get_spin(uint32_t i) const { 
      Spin3D s = lattice.spins[i];
      if (cluster_update) {
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

    void cluster_mutation();
    void metropolis_mutation();

    virtual double energy() const override;
    virtual double energy_change() override;
    virtual void generate_mutation() override;
    virtual void accept_mutation() override;
    virtual void reject_mutation() override;

    virtual double onsite_func(const Spin3D& S) const=0;
    virtual double onsite_energy(uint32_t i) const;
    virtual double bond_energy(uint32_t i) const;

    void add_magnetization_samples(dataframe::SampleMap &samples) const;
    void add_helicity_samples(dataframe::SampleMap &samples) const;
    void add_spin_samples(dataframe::SampleMap& samples) const;

    void add_intensityx_samples(dataframe::SampleMap& samples) const;
    void add_intensityy_samples(dataframe::SampleMap& samples) const;
    void add_intensityz_samples(dataframe::SampleMap& samples) const;

    virtual dataframe::SampleMap take_samples() override;

    LatticeGraph to_graph() const { 
      return lattice.to_graph();
    }

  protected:
    Lattice<Spin3D> lattice;

    struct Spin3DMutation {
      uint32_t i;
      Spin3D dS;
    };

    bool cluster_update;

    // Mutation being considered is stored as an attribute of the model
    Spin3DMutation mut;

    static constexpr double alpha = 0.01;

  private:
    bool sample_energy;
    bool sample_magnetization;
    bool sample_helicity;
    bool sample_spins;

    bool sample_intensityx;
    bool sample_intensityy;
    bool sample_intensityz;
    double max_L;
    double min_L;
    uint32_t intensity_resolution;

    uint64_t nsteps;
    uint64_t accepted;
    double acceptance;
    double sigma;

    std::unordered_set<int> s;
    Eigen::Matrix3d s0;

    // Internal normally distributed random number generator
    GaussianDist dist;
};

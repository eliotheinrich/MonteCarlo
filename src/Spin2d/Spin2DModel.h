#pragma once

#include "MonteCarlo.h"

#include <Eigen/Dense>

#include <functional>
#include <vector>
#include <unordered_set>
#include <stack>

class Spin2DModel : public MonteCarloSimulator {
  // Generic 2d spin model in up to 3d lattice
  public:
    uint32_t V;

    Spin2DModel(dataframe::Params &params, uint32_t num_threads);
    virtual ~Spin2DModel()=default;

    void init(uint32_t sl, uint32_t N1, uint32_t N2, uint32_t N3);

    virtual uint64_t system_size() const override {
      if (cluster_update) {
        return 1;
      }

      return V;
    }

    inline uint32_t flat_idx(int n1, int n2, int n3, int s) const {
      return n1 + N1*(n2 + N2*(n3 + N3*s));
    }

    inline Eigen::Vector4i tensor_idx(int i) const {
      int n1 = i % N1;
      i = i / N1;
      int n2 = i % N2;
      i = i / N2;
      int n3 = i % N3;
      i = i / N3;
      int s = i % sl;
      Eigen::Vector4i v; v << n1, n2, n3, s;
      return v;
    }

    void randomize_spins();

    Eigen::Vector2d get_spin(uint32_t i) const {
      return cluster_update ? s0.transpose()*spins[i] : spins[i];
    }

    void add_bond(
        int d1, 
        int d2, 
        int d3, 
        int ds, 
        Eigen::Vector3d v, 
        std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)> bondfunc
    ) {
        add_bond(d1, d2, d3, ds, v, bondfunc, [](uint32_t, uint32_t) { return true; });
    }

    void add_bond(
        int d1, 
        int d2, 
        int d3, 
        int ds, 
        Eigen::Vector3d v, 
        std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)> bondfunc, 
        std::function<bool(uint32_t, uint32_t)> bond_filter
    );

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

    void add_magnetization_samples(dataframe::data_t &samples) const;
    void add_helicity_samples(dataframe::data_t &samples) const;
    virtual dataframe::data_t take_samples() override;

    void save_spins(const std::string& filename);

  protected:
    // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
    // dS must conserve the norm of S[n1,n2,n3,s]
    struct Spin2DMutation {
      uint32_t i;
      Eigen::Vector2d dS;
    };

    struct Spin2DBond {
      int d1;
      int d2;
      int d3;
      int ds;
      Eigen::Vector3d v;
      std::function<double(Eigen::Vector2d, Eigen::Vector2d)> bondfunc;
    };

    static constexpr double alpha = 0.01;
    std::vector<Eigen::Matrix2d> R1s;
    std::vector<Eigen::Matrix2d> R2s;
    std::vector<Eigen::Matrix2d> R3s;

    bool cluster_update;

    std::vector<std::vector<Bond>> neighbors;
    std::vector<Spin2DBond> bonds;

    // Mutation being considered is stored as an attribute of the model
    Spin2DMutation mut;

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

    double acceptance;
    double sigma;
    std::vector<Eigen::Vector2d> spins;

    std::vector<std::function<bool(uint32_t, uint32_t)>> bond_filters;

    std::unordered_set<uint32_t> s;
    Eigen::Matrix2d s0;
};

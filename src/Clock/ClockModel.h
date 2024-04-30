#pragma once

#include "MonteCarlo.h"

#include <Eigen/Dense>

#include <vector>
#include <unordered_set>
#include <stack>
#include <math.h>
#include <fstream>

template <uint32_t q>
class ClockModel : public MonteCarloSimulator {
  // Generic 3D Ising model
  public:
    uint32_t N1;
    uint32_t N2;
    uint32_t N3;
    uint64_t V;

    std::vector<int> spins;
    std::vector<std::vector<uint32_t>> neighbors;

    std::unordered_set<uint32_t> s;

    // Mutation being considered is stored as an attribute of the model
    int mut_mode;

    ClockModel(dataframe::Params& params, uint32_t num_threads) : MonteCarloSimulator(params, num_threads) {}
    ClockModel()=default;
    virtual ~ClockModel()=default;

    void init(uint32_t N1, uint32_t N2, uint32_t N3) {
      this->N1 = N1;
      this->N2 = N2;
      this->N3 = N3;
      this->V = N1*N2*N3;
    }

    virtual void init_state() override {
      spins = std::vector<int>(V);
      neighbors = std::vector<std::vector<uint32_t>>(V, std::vector<uint32_t>(0));
      randomize_spins();

      mut.i = 0;
      mut_mode = 0;
    }

    virtual uint64_t system_size() const override {
      return V;
    }

    inline const uint32_t flat_idx(uint32_t n1, uint32_t n2, uint32_t n3) {
      return n1 + N1*(n2 + N2*n3);
    }

    inline const Eigen::Vector3i tensor_idx(uint32_t i) {
      uint32_t n1 = i % N1;
      i = i / N1;
      uint32_t n2 = i % N2;
      i = i / N2;
      uint32_t n3 = i % N3;
      Eigen::Vector3i v; v << n1, n2, n3;
      return v;
    }

    void randomize_spins() {
      for (uint32_t i = 0; i < V; i++) {
        // For each site, initialize spin randomly
        spins[i] = rand() % q;
      }
    }

    void add_bond(
      int d1, 
      int d2, 
      int d3, 
      Eigen::Vector3d v, 
      std::function<double(uint32_t, uint32_t)> bondfunc
    ) {
      ClockBond b{d1, d2, d3, v, bondfunc};
      this->bonds.push_back(b);
      for (uint32_t n1 = 0; n1 < N1; n1++) {
        for (uint32_t n2 = 0; n2 < N2; n2++) {
          for (uint32_t n3 = 0; n3 < N3; n3++) {
            uint32_t i = flat_idx(n1, n2, n3);
            uint32_t j = flat_idx(mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3));
            neighbors[i].push_back(j);
          }
        }
      }
    }

    virtual double onsite_energy(uint32_t i) const=0;

    virtual double bond_energy(uint32_t i) const {
      double E = 0.;
      for (uint32_t n = 0; n < bonds.size(); n++) {
        uint32_t j = neighbors[i][n];
        E += 0.5*bonds[n].bondfunc(spins[i], spins[j]);
      }

      return E;
    }

    virtual double energy() const override {
      double E = 0;

      for (uint32_t i = 0; i < V; i++) {
        E += onsite_energy(i);
        E += bond_energy(i);
      }
      return E;
    }

    inline double get_magnetization() const {
      double x = 0; 
      double y = 0;
      for (uint32_t i = 0; i < V; i++) {
        x += std::cos(2*PI*spins[i]/q);
        y += std::sin(2*PI*spins[i]/q);
      }

      return std::sqrt(x*x + y*y)/(N1*N2*N3);
    }

    void metropolis_mutation() {
      mut.dq = rand() % 3 - 1;
    }

    void cluster_mutation() {
      s.clear();

      std::stack<int> c;
      int m = rand() % V;
      c.push(m);

      int p = rand() % q;

      while (!c.empty()) {
        m = c.top();
        c.pop();

        if (!s.count(m)) {
          s.insert(m);

          uint32_t s_new = mod(2*spins[m] - p, q);

          for (uint32_t n = 0; n < neighbors[m].size(); n++) {
            uint32_t j = neighbors[m][n];

            if (!s.count(j)) {
              double dE = bonds[n].bondfunc(spins[j], s_new) - bonds[n].bondfunc(spins[j], spins[m]);
              if (randf() < 1. - std::exp(-dE/temperature)) {
                c.push(j);
              }
            }
          }
          spins[m] = s_new;
        }
      }
    }

    virtual void generate_mutation() override {
      mut.i++;
      if (mut.i == V) {
        mut.i = 0;
        mut_mode++;
      }

      if (mut_mode < 3) {
        metropolis_mutation();
      } else {
        cluster_mutation();
        mut_mode = 0;
      }
    }

    virtual double energy_change() override {
      if (mut_mode == 0) { 
        return -1.;
      }
      else {
        double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
        spins[mut.i] = mod(spins[mut.i] + mut.dq, q);
        double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
        return E2 - E1;
      }
    }

    virtual void accept_mutation() override {
      return;
    }

    virtual void reject_mutation() override {
      spins[mut.i] = mod(spins[mut.i] - mut.dq, q);
    }


    // Saves current spin configuration
    void save_spins(const std::string& filename) {
      std::ofstream output_file(filename);

      output_file << N1 << "\t" << N2 << "\t" << N3 << "\n";
      for (uint32_t i = 0; i < V; i++) {
        output_file << cos(2*PI*spins[i]/q) << "\t" << sin(2*PI*spins[i]/q);
        if (i < V - 1) { 
          output_file << "\t"; 
        }
      }
      output_file.close();
    }

    virtual dataframe::data_t take_samples() const {
      dataframe::data_t samples;
      samples.emplace("energy", energy());
      samples.emplace("magnetization", get_magnetization());
      return samples;
    }

  private:
  // Must be supplied with number of sublattices
    struct ClockMutation {
      int dq;
      uint32_t i;
    };

    struct ClockBond {
      int d1;
      int d2;
      int d3;
      Eigen::Vector3d v;
      std::function<double(uint32_t, uint32_t)> bondfunc;
    };    

  protected:
    std::vector<ClockBond> bonds;
    ClockMutation mut;
};


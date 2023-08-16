#ifndef SPIN3D_MC_H
#define SPIN3D_MC_H

#include <vector>
#include <unordered_set>
#include <math.h>
#include <random>
#include <Eigen/Dense>
#include "MonteCarlo.h"

class GaussianDist {
    private:
        std::minstd_rand rd;
        std::default_random_engine gen;
        std::normal_distribution<> dist;

    public:
        GaussianDist() {}
        GaussianDist(float mean, float std);

        float sample();

};

class Spin3DModel : virtual public MCModel {
    // Generic 3D Heisenberg model
    private:
        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct Spin3DMutation {
            int i;
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

        bool sample_energy;
        bool sample_magnetization;
        bool sample_helicity;

        int sl;
        int N1;
        int N2;
        int N3;

        BoundaryCondition bcx;
        BoundaryCondition bcy;
        BoundaryCondition bcz;

        ull nsteps;
        ull accepted;
        float acceptance;
        double sigma;
        std::vector<Eigen::Vector3d> spins;

        // Need to store bond filters passed to add_bond so that neighbor matrix can 
        // be later filled out bit init
        std::vector<std::function<bool(uint, uint)>> bond_filters;

        std::unordered_set<int> s;
        Eigen::Matrix3d s0;

        // Internal normally distributed random number generator
        GaussianDist dist;

    protected:
        bool cluster_update;

        // Mutation being considered is stored as an attribute of the model
        Spin3DMutation mut;

        std::vector<std::vector<Bond>> neighbors;
        std::vector<Spin3DBond> bonds;

        static constexpr double alpha = 0.01;
        std::vector<Eigen::Matrix3d> R1s;
        std::vector<Eigen::Matrix3d> R2s;
        std::vector<Eigen::Matrix3d> R3s;


    public:
        ull V;

        Spin3DModel(Params &params);
        virtual ~Spin3DModel() {}

        void init_params(int sl, int N1, int N2, int N3);
        virtual void init() override;

        virtual ull system_size() const override {
            if (cluster_update)
                return 1;
            
            return V;
        }

        void set_spin(int i, Eigen::Vector3d S) { spins[i] = S; }
        Eigen::Vector3d get_spin(int i) const { return cluster_update ? s0.transpose()*spins[i] : spins[i]; }
        void randomize_spins();

        int flat_idx(int n1, int n2, int n3, int s) const {
            return n1 + N1*(n2 + N2*(n3 + N3*s));
        }

        Eigen::Vector4i tensor_idx(int i) const {
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

        void add_bond(int d1, int d2, int d3, int ds, 
                      Eigen::Vector3d v, 
                      std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc) {
            add_bond(d1, d2, d3, ds, v, bondfunc, [](uint, uint){ return true; });
        }
        void add_bond(int d1, int d2, int d3, int ds, 
                      Eigen::Vector3d v, 
                      std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc,
                      std::function<bool(uint, uint)> bond_filter);

        static std::vector<double> twist_terms(std::vector<double> dE);
        std::vector<double> twist_derivatives(int i) const;
        std::vector<double> twist_derivatives() const;

        Eigen::Vector3d get_magnetization() const;

        std::vector<double> correlation_function(int i, int a, int b) const;
        std::vector<double> full_correlation_function(int i) const;
        double skyrmion_density(int i) const;

        std::vector<double> skyrmion_correlation_function(int i) const;

        void cluster_mutation();
        void metropolis_mutation();

        virtual double energy() const override;
        virtual double energy_change() override;
        virtual void generate_mutation() override;
        virtual void accept_mutation() override;
        virtual void reject_mutation() override;

        virtual double onsite_func(const Eigen::Vector3d& S) const = 0;
        virtual double onsite_energy(int i) const;
        virtual double bond_energy(int i) const;


        // Saves current spin configuration
        void save_spins(std::string filename);
        bool load_spins(std::string filename);

        void add_magnetization_samples(data_t &samples) const;
        void add_helicity_samples(data_t &samples) const;
        virtual data_t take_samples() override;
};

#endif

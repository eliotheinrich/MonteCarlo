#ifndef SPINMC_H
#define SPINMC_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <stack>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <fstream>
#include <Eigen/Dense>
#include "MonteCarlo.h"
#include "SpinModel.h"
#include "Utility.h"

class SpinModel : virtual public MCModel {
    // Generic 3D Heisenberg model
    private:
        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct SpinMutation {
            int i;
            Eigen::Vector3d dS;
        };

        struct HeisBond {
            int d1;
            int d2;
            int d3;
            int ds;
            Eigen::Vector3d v;
            std::function<double(const Eigen::Vector3d&, const Eigen::Vector3d&)> bondfunc;
        };    


    public:
        int sl;
        int N1;
        int N2;
        int N3;
        int V;

        float acceptance;
        double sigma;
        std::vector<Eigen::Vector3d> spins;
        std::vector<std::vector<int>> neighbors;
        std::vector<HeisBond> bonds;

        static constexpr double alpha = 0.01;
        std::vector<Eigen::Matrix3d> R1s;
        std::vector<Eigen::Matrix3d> R2s;
        std::vector<Eigen::Matrix3d> R3s;

        std::unordered_set<int> s;
        Eigen::Matrix3d s0;

        bool tracking;
        std::vector<double> q;

        std::minstd_rand r;

        // Mutation being considered is stored as an attribute of the model
        SpinMutation mut;

        // Internal normally distributed random number generator
        GaussianDist *dist;

        SpinModel();

        SpinModel(int sl, int N1, int N2, int N3);

        virtual std::vector<double> tracking_func(int i);
        virtual std::vector<double> init_func();
        void start_tracking();

        void set_spin(int i, Eigen::Vector3d S);
        void randomize_spins();

        inline int flat_idx(int n1, int n2, int n3, int s) const {
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

        void add_bond(int d1, int d2, int d3, int ds, 
                      Eigen::Vector3d v, 
                      std::function<double(Eigen::Vector3d, Eigen::Vector3d)> bondfunc);

        static std::vector<double> twist_terms(std::vector<double> dE);
        std::vector<double> twist_derivatives(int i);
        const std::vector<double> twist_derivatives();

        Eigen::Vector3d get_magnetization() const;

        const std::vector<double> correlation_function(int i, int a, int b);
        const std::vector<double> full_correlation_function(int i);
        const double skyrmion_density(int i);

        const std::vector<double> skyrmion_correlation_function(int i);

        void cluster_update();
        void generate_mutation();
        void metropolis_mutation();

        void accept_mutation();
        void reject_mutation();

        virtual double onsite_func(const Eigen::Vector3d& S) const = 0;
        virtual double onsite_energy(int i) const;
        virtual double bond_energy(int i) const;

        double energy() const;
        double energy_change();

        // Saves current spin configuration
        void save_spins(std::string filename);
        bool load_spins(std::string filename);

        virtual std::map<std::string, double> get_double_params() const;
        virtual std::map<std::string, int> get_int_params() const;
};

#endif

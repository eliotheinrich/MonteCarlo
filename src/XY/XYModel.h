#ifndef XYMC_H
#define XYMC_H

#include <iostream>
#include <functional>
#include <vector>
#include <unordered_set>
#include <stdlib.h>
#include <math.h>
#include <stack>
#include <random>
#include <Eigen/Dense>
#include <fstream>
#include "MonteCarlo.h"
#include "Utility.h"
    
class XYModel : virtual public MCModel {
    // Generic 3D XY model
    private:
        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct XYMutation {
            int i;
            Eigen::Vector2d dS;
        };

        struct XYBond {
            int d1;
            int d2;
            int d3;
            int ds;
            Eigen::Vector3d v;
            std::function<float(Eigen::Vector2d, Eigen::Vector2d)> bondfunc;
        };



    public:
        int sl;
        int N1;
        int N2;
        int N3;
        int V;

        float acceptance;
        float sigma;
        std::vector<Eigen::Vector2d> spins;
        std::vector<std::vector<int>> neighbors;
        std::vector<XYBond> bonds;

        static constexpr double alpha = 0.01;
        std::vector<Eigen::Matrix2d> R1s;
        std::vector<Eigen::Matrix2d> R2s;
        std::vector<Eigen::Matrix2d> R3s;


        std::unordered_set<int> s;
        Eigen::Matrix2d s0;

        std::minstd_rand r;

        // Mutation being considered is stored as an attribute of the model
        XYMutation mut;

        XYModel() {}

        XYModel(int sl, int N1, int N2, int N3);

        inline const int flat_idx(int n1, int n2, int n3, int s) {
            return n1 + N1*(n2 + N2*(n3 + N3*s));
        }

        inline const Eigen::Vector4i tensor_idx(int i) {
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

        void add_bond(int d1, int d2, int d3, int ds, Eigen::Vector3d v, std::function<float(Eigen::Vector2d, Eigen::Vector2d)> bondfunc);

        std::vector<double> twist_stiffness();

        Eigen::Vector2d get_magnetization();
        
        virtual double onsite_energy(int i) const;
        virtual double bond_energy(int i) const;
        virtual double energy() const;

        void generate_mutation();
        void metropolis_mutation();
        void cluster_update();

        void accept_mutation();
        void reject_mutation();
        virtual double energy_change();

        virtual double onsite_func(const Eigen::Vector2d& S) const = 0;


        void save_spins(std::string filename);
};

#endif

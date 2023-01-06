#ifndef CLOCKMC_H
#define CLOCKMC_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <stack>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <Eigen/Dense>
#include <fstream>
#include "MonteCarlo.h"
#include "Utility.cpp"

template <int q>
class ClockModel : virtual public MCModel {
    // Generic 3D Ising model
    // Must be supplied with number of sublattices
    private:
        struct ClockMutation {
            int dq;
            int i;
        };

        struct ClockBond {
            int d1;
            int d2;
            int d3;
            Eigen::Vector3d v;
            std::function<float(int, int)> bondfunc;
        };    


    public:
        int N1;
        int N2;
        int N3;
        int V;

        std::vector<int> spins;
        std::vector<ClockBond> bonds;
        std::vector<std::vector<int>> neighbors;

        std::unordered_set<int> s;

        std::minstd_rand r;

        // Mutation being considered is stored as an attribute of the model
        ClockMutation mut;
        int mut_mode;

        ClockModel();
        ClockModel(int N1, int N2, int N3);

        inline const int flat_idx(int n1, int n2, int n3);

        inline const Eigen::Vector3i tensor_idx(int i);

        void randomize_spins();

        void add_bond(int d1, int d2, int d3, Eigen::Vector3d v, std::function<float(int, int)> bondfunc);

        virtual const float onsite_energy(int i)=0;
        virtual const float bond_energy(int i);
        const float energy();
        inline float get_magnetization();

        void metropolis_mutation();
        void cluster_update();

        void generate_mutation();
        const float energy_change();
        void accept_mutation();
        void reject_mutation();


        // Saves current spin configuration
        void save_spins(std::string filename);
};

#include "ClockModel.cpp"

#endif

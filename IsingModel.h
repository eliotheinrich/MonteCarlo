#ifndef ISINGMC_H
#define ISINGMC_H

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include "MonteCarlo.h"
#include "Utility.cpp"

class IsingModel : virtual public MCModel {
    // Generic 3D Ising model
    // Must be supplied with number of sublattices
    private:
        struct IsingMutation {
            int i;
        };


    public:
        int N1;
        int N2;
        int N3;
        int V;

        float acceptance;
        std::vector<float> spins;

        // Mutation being considered is stored as an attribute of the model
        IsingMutation mut;

        IsingModel() {}

        IsingModel(int N1, int N2, int N3);

        inline int flat_idx(int n1, int n2, int n3) const;

        inline Eigen::Vector3i tensor_idx(int i) const;

        void randomize_spins();

        inline double get_magnetization() const;
        virtual double energy() const;
        virtual double onsite_energy(int i) const = 0;
        virtual double bond_energy(int i) const = 0;

        virtual void generate_mutation();
        virtual double energy_change();
        virtual void accept_mutation();
        virtual void reject_mutation();

        void save_spins(std::string filename);
};

#endif

#ifndef ISINGMC_H
#define ISINGMC_H

#include <vector>
#include <Eigen/Dense>
#include "MonteCarlo.h"

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

        void init_params(int N1, int N2, int N3);

        virtual void init();

        inline int flat_idx(int n1, int n2, int n3) const {
            return n1 + N1*(n2 + N2*n3);
        }

        inline Eigen::Vector3i tensor_idx(int i) const {
            int n1 = i % N1;
            i = i / N1;
            int n2 = i % N2;
            i = i / N2;
            int n3 = i % N3;
            Eigen::Vector3i v; v << n1, n2, n3;
            return v;
        }

        void randomize_spins();

        double get_magnetization() const;
        virtual double onsite_energy(int i) const = 0;
        virtual double bond_energy(int i) const = 0;

        virtual void generate_mutation();
        virtual void accept_mutation();
        virtual void reject_mutation();
        
        virtual double energy() const;
        virtual double energy_change();

        virtual std::map<std::string, Sample> take_samples() const;

        void save_spins(std::string filename);
};

#endif

#ifndef SPINMC_
#define SPINMC_

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <regex>
#include "MonteCarlo.cpp"
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

        IsingModel(int N1, int N2, int N3) {
            this->N1 = N1;
            this->N2 = N2;
            this->N3 = N3;
            this->V = N1*N2*N3;

            this->spins = std::vector<float>(V);

            this->randomize_spins();

            this->acceptance = 0.5;
        }

        inline const int flat_idx(int n1, int n2, int n3) {
            return n1 + N1*(n2 + N2*n3);
        }

        inline const Eigen::Vector3i tensor_idx(int i) {
            int n1 = i % N1;
            i = i / N1;
            int n2 = i % N2;
            i = i / N2;
            int n3 = i % N3;
            Eigen::Vector3i v; v << n1, n2, n3;
            return v;
        }

        void randomize_spins() {
            for (int i = 0; i < V; i++) {
                // For each site, initialize spin randomly
                if (rand() % 2) {
                    this->spins[i] = 1.;
                } else {
                    this->spins[i] = -1.;
                }
            }
        }

        inline float get_magnetization() {
            float M = 0;
            for (int i = 0; i < V; i++) {
                M += this->spins[i];
            }
            
            return M/(N1*N2*N3);
        }

        void generate_mutation() {
            // Randomly select a site to mutate
            mut.i = std::rand() % V;
        }

        void accept_mutation() {
            return;
        }

        void reject_mutation() {
            this->spins[mut.i] = -this->spins[mut.i];
        }

        virtual const float onsite_energy(int i)=0;

        virtual const float bond_energy(int i)=0;

        const float energy() {
            float E = 0;

            for (int i = 0; i < V; i++) {
                E += onsite_energy(i);
                E += bond_energy(i);
            }
            return E;
        }

        const float energy_change() {
            float E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
            this->spins[mut.i] = -this->spins[mut.i];
            float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

            return E2 - E1;
        }

        // Saves current spin configuration
        void save_spins(std::string filename) {
            std::ofstream output_file;
            output_file.open(filename);
            output_file << N1 << "\t" << N2 << "\t" << N3 << "\n";
            for (int i = 0; i < V; i++) {
                output_file << spins[i];
                if (i < V-1) { output_file << "\t"; }
            }
            output_file.close();
        }        
};

#endif

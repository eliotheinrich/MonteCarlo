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

using namespace std;
using namespace Eigen;


class IsingModel : virtual public MCModel {
    // Generic 3D Ising model
    // Must be supplied with number of sublattices
    private:
        struct IsingMutation {
            int n1;
            int n2;
            int n3;
        };


    public:
        int N1;
        int N2;
        int N3;
        int V;

        float acceptance;
        vector<vector<vector<float>>> spins;

        // Mutation being considered is stored as an attribute of the model
        LatticeIterator* iter;
        IsingMutation mut;

        IsingModel() {}

        IsingModel(int N1, int N2, int N3) {
            this->N1 = N1;
            this->N2 = N2;
            this->N3 = N3;
            this->V = N1*N2*N3;

            this->spins = vector<vector<vector<float>>>(this->N1, vector<vector<float>>(this->N2, vector<float>(this->N3)));

            this->iter = new LatticeIterator(N1, N2, N3, 1);

            this->randomize_spins();

            this->acceptance = 0.5;
        }

        void randomize_spins() {
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        // For each site, initialize spin randomly
                        if (rand() % 2) {
                            this->spins[n1][n2][n3] = 1.;
                        } else {
                            this->spins[n1][n2][n3] = -1.;
                        }
                    }
                }
            }
        }

        inline float get_magnetization() {
            float M = 0;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        M += this->spins[n1][n2][n3];
                    }
                }
            }
            
            return M/(N1*N2*N3);
        }

        void generate_mutation() {
            // Randomly select a site to mutate
            int n1 = iter->n1;
            int n2 = iter->n2;
            int n3 = iter->n3;
            iter->next();

            this->mut.n1 = n1;
            this->mut.n2 = n2;
            this->mut.n3 = n3;
        }

        void accept_mutation() {
            return;
        }

        void reject_mutation() {
            this->spins[mut.n1][mut.n2][mut.n3] = -this->spins[mut.n1][mut.n2][mut.n3];
        }

        virtual const float onsite_energy(int n1, int n2, int n3)=0;

        virtual const float bond_energy(int n1, int n2, int n3)=0;

        const float energy() {
            float E = 0;

            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        E += onsite_energy(n1, n2, n3);
                        E += bond_energy(n1, n2, n3);
                    }
                }
            }
            return E;
        }

        const float energy_change() {
            float E1 = onsite_energy(mut.n1, mut.n2, mut.n3) + 2*bond_energy(mut.n1, mut.n2, mut.n3);
            this->spins[mut.n1][mut.n2][mut.n3] = -this->spins[mut.n1][mut.n2][mut.n3];
            float E2 = onsite_energy(mut.n1, mut.n2, mut.n3) + 2*bond_energy(mut.n1, mut.n2, mut.n3);

            return E2 - E1;
        }

        // Saves current spin configuration
        void save_spins(string filename) {
            ofstream output_file;
            output_file.open(filename);
            output_file << N1 << "\t" << N2 << "\t" << N3 << endl;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        output_file << spins[n1][n2][n3] << ",";
                    }
                    output_file << ";";
                }
                output_file << endl;
            }
            output_file.close();
        }        
};

#endif

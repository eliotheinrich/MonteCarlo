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
            Vector3f v;
            function<float(int, int)> bondfunc;
        };    


    public:
        int N1;
        int N2;
        int N3;
        int V;

        vector<int> spins;
        vector<ClockBond> bonds;
        vector<vector<int>> neighbors;

        minstd_rand r;

        // Mutation being considered is stored as an attribute of the model
        ClockMutation mut;

        ClockModel() {}

        ClockModel(int N1, int N2, int N3) {
            this->N1 = N1;
            this->N2 = N2;
            this->N3 = N3;
            this->V = N1*N2*N3;

            this->spins = vector<int>(V);
            this->neighbors = vector<vector<int>>(V, vector<int>(0));

            this->r.seed(rand());
            this->randomize_spins();

            this->mut.i = 0;
        }

        inline const int flat_idx(int n1, int n2, int n3) {
            return n1 + N1*(n2 + N2*n3);
        }

        inline const Vector3i tensor_idx(int i) {
            int n1 = i % N1;
            i = i / N1;
            int n2 = i % N2;
            i = i / N2;
            int n3 = i % N3;
            Vector3i v; v << n1, n2, n3;
            return v;
        }

        void randomize_spins() {
            for (int i = 0; i < V; i++) {
                // For each site, initialize spin randomly
                spins[i] = r() % q;
            }
        }

        void add_bond(int d1, int d2, int d3, Vector3f v, function<float(int, int)> bondfunc) {
            ClockBond b{d1, d2, d3, v, bondfunc};
            this->bonds.push_back(b);
            int i; int j;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        i = flat_idx(n1, n2, n3);
                        j = flat_idx(mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3));
                        neighbors[i].push_back(j);
                    }
                }
            }
        }

        inline float get_magnetization() {
            float M = 0;
            float x = 0; float y = 0;
            for (int i = 0; i < V; i++) {
                x += cos(2*PI*spins[i]/q);
                y += sin(2*PI*spins[i]/q);
            }
            
            return sqrt(x*x + y*y)/(N1*N2*N3);
        }

        void generate_mutation() {
            // Randomly select a site to mutate
            mut.i++;
            if (mut.i == V) {
                mut.i = 0;
            }
            mut.dq = r() % 3 - 1;
        }

        void accept_mutation() {
            return;
        }

        void reject_mutation() {
            this->spins[mut.i] = mod(this->spins[mut.i] - mut.dq, q);
        }

        virtual const float onsite_energy(int i)=0;

        virtual const float bond_energy(int i) {
            float E = 0.;
            int j;
            for (int n = 0; n < bonds.size(); n++) {
                j = neighbors[i][n];
                E += 0.5*bonds[n].bondfunc(spins[i], spins[j]);
            }

            return E;
        }

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
            this->spins[mut.i] = mod(this->spins[mut.i] + mut.dq, q);
            float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

            return E2 - E1;
        }

        // Saves current spin configuration
        void save_spins(string filename) {
            ofstream output_file;
            output_file.open(filename);
            output_file << N1 << "\t" << N2 << "\t" << N3 << "\n";
            for (int i = 0; i < V; i++) {
                output_file << cos(2*PI*spins[i]/q) << "\t" << sin(2*PI*spins[i]/q);
                if (i < V-1) { output_file << "\t"; }
            }
            output_file.close();
        }        
};

#endif

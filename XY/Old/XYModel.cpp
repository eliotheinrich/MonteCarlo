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


//template<unsigned int sublattices>
class XYModel : virtual public MCModel {
    // Generic 3D Heisenberg model
    // Must be supplied with number of sublattices
    private:
        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct XYMutation {
            int n1;
            int n2;
            int n3;
            int s;
            float dp;
        };


    public:
        int sl;
        int N1;
        int N2;
        int N3;
        float acceptance;
        float sigma;
        vector<vector<vector<vector<float>>>> spins;

        // Mutation being considered is stored as an attribute of the model
        XYMutation mut;

        XYModel() {}

        XYModel(int sl, int N1, int N2 = -1, int N3 = -1) {
            this->sl = sl;
            this->N1 = N1;
            if (N2 == -1) { this->N2 = N1; } else { this->N2 = N2; }
            if (N3 == -1) { this->N3 = N1; } else { this->N3 = N3; }
            this->spins = vector<vector<vector<vector<float>>>>(this->N1,
                                 vector<vector<vector<float>>>(this->N2,
                                        vector<vector<float>>(this->N3,
                                               vector<float>(sl)))); 

            this->randomize_spins();

            this->acceptance = 0.5;
            this->sigma = 0.25;
        }

        void randomize_spins() {
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            // For each site, initialize spin randomly
                            spins[n1][n2][n3][s] = 2*PI*float(rand())/float(RAND_MAX)
                        }
                    }
                }
            }
        }

        inline Vector2f get_magnetization(int s) {
            Vector2f M = Vector2f::Constant(0);
            Vector2f tmp;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        tmp << cos(spins[n1][n2][n3][s]), sin(spins[n1][n2][n3][s]);
                        M += tmp;
                    }
                }
            }
            
            return M/(N1*N2*N3);
        }

        inline Vector2f get_magnetization() {
            Vector2f M = Vector2f::Constant(0);
            for (int s = 0; s < sl; s++) {
                M += get_magnetization(s);
            }
            return M/(sl);
        }

        void generate_mutation() {
            // TODO
            // Find a better way to update cone width
            if (acceptance > 0.5) {
                sigma = min(2., 1.01*sigma);

            } else {
                sigma = max(0.05, 0.99*sigma);
            }

            // Randomly select a site to mutate
            int rand_n1 = rand() % N1;
            int rand_n2 = rand() % N2;
            int rand_n3;
            if (this->N3 == 1) { rand_n3 = 0; } else { rand_n3 = rand() % N3; }

            int rand_s;
            if (sl == 1) { rand_s = 0; } else { rand_s = rand() % sl; }

            float dp = sigma*2*PI*float(rand())/float(RAND_MAX);

            // Store mutation for consideration
            this->mut.n1 = rand_n1;
            this->mut.n2 = rand_n2;
            this->mut.n3 = rand_n3;
            this->mut.s = rand_s;
            this->mut.dp = dp;
        }

        void accept_mutation() {
            return;
        }

        void reject_mutation() {
            this->spins[mut.n1][mut.n2][mut.n3][mut.s] -= mut.dp;
        }

        const float onsite_energy(int n1, int n2, int n3, int s)=0;

        const float bond_energy(int n1, int n2, int n3, int s)=0;

        const float energy() {
            float E = 0;

            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            E += onsite_energy(n1, n2, n3, s);
                            E += bond_energy(n1, n2, n3, s);
                        }
                    }
                }
            }
            return E;
        }

        const float energy_change() {
            float E1 = onsite_energy(mut.n1, mut.n2, mut.n3, mut.s) + 2*bond_energy(mut.n1, mut.n2, mut.n3, mut.s);
            this->spins[mut.n1][mut.n2][mut.n3][mut.s] += mut.dp;
            float E2 = onsite_energy(mut.n1, mut.n2, mut.n3, mut.s) + 2*bond_energy(mut.n1, mut.n2, mut.n3, mut.s);

            return E2 - E1;
        }

        // Saves current spin configuration
        void save_spins(string filename) {
            ofstream output_file;
            output_file.open(filename);
            output_file << sl << "\t" << N1 << "\t" << N2 << "\t" << N3 << endl;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s<sl; s++) {
                            output_file << "[" << spins[n1][n2][n3][s] << "]";
                            if (s != sl - 1) {
                                output_file << "\t";
                            }
                        }
                        output_file << ",";
                    }
                    output_file << ";";
                }
                output_file << endl;
            }
            output_file.close();
        }        
};

//template<int sublattices>
class MagnetizationLogItem : public LogItem {
    // Stores total magnetization, acceptance rate, and energy
    public:
        Vector2f magnetization;
        float acceptance;
        float energy;

        MagnetizationLogItem() {
            magnetization = Vector2f::Constant(0);
            acceptance = 0.;
            energy = 0.;
        }

        MagnetizationLogItem(MonteCarlo *m, XYModel *model) {
            magnetization = model->get_magnetization();
            acceptance = model->acceptance;
            energy = m->energy;
        }

        friend ostream& operator<<(ostream& os, const MagnetizationLogItem& logitem) {
            os << logitem.magnetization[0] << " " << logitem.magnetization[1] << " " << logitem.magnetization[2] << "\t"
               << logitem.acceptance << "\t"
               << logitem.energy;
            return os;
        }
};


#endif

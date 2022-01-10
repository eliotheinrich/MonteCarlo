#ifndef SPINMC_
#define SPINMC_

#include <iostream>
#include <functional>
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

struct Bond {
    int d1;
    int d2;
    int d3;
    int ds;
    Vector3f v;
    function<float(Vector2f, Vector2f)> bondfunc;
};
    
class XYModel : virtual public MCModel {
    // Generic 3D XY model
    private:
        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct XYMutation {
            int i;
            Vector2f dS;
        };


    public:
        int sl;
        int N1;
        int N2;
        int N3;
        int V;

        float acceptance;
        float sigma;
        vector<Vector2f> spins;
        vector<vector<int>> neighbors;
        vector<Bond> bonds;

        minstd_rand r;

        int mut_counter;
        int mut_mode;

        // Mutation being considered is stored as an attribute of the model
        XYMutation mut;

        XYModel() {}

        XYModel(int sl, int N1, int N2 = -1, int N3 = -1) {
            this->sl = sl;
            this->N1 = N1;
            if (N2 == -1) { this->N2 = N1; } else { this->N2 = N2; }
            if (N3 == -1) { this->N3 = N1; } else { this->N3 = N3; }
            this->V = N1*N2*N3*sl;

            this->spins = vector<Vector2f>(V);
            this->neighbors = vector<vector<int>>(V, vector<int>(0));

            this->randomize_spins();

            this->acceptance = 0.5;
            this->sigma = 0.25;

            this->r.seed(rand());

            this->mut_counter = 0;
            this->mut_mode = 0;
        }

        void randomize_spins() {
            float p;
            Vector2f v;

            for (int i = 0; i < V; i++) {
                p = 2*PI*float(r())/float(RAND_MAX);
                v << cos(p), sin(p);
                set_spins(i, v);
            }
        }

        void add_bond(Bond b) {
            Vector4i v;
            this->bonds.push_back(b);
            int idx; int neighbor_idx;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            idx = flat_idx(n1, n2, n3, s);
                            neighbor_idx = flat_idx(mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3), mod(s + b.ds, sl));
                            neighbors[i].push_back(neighbor_idx);
                        }
                    }
                }
            }
        }

        inline vector<double> twist_stiffness() {
            // Returns the first and second derivative in response to a phase twist
            float alpha = 0.01;
            float f;
            Matrix2f R1; 
            Matrix2f R2;

            double E0 = 0.;
            double E1 = 0.;
            double E2 = 0.;
            double Em1 = 0.;
            double Em2 = 0.;


            Vector2f S1;
            Vector2f S2;
            int j;
            for (int i = 0; i < V; i++) {
                for (int n = 0; n < bonds.size(); n++) {
                    f = bonds[0].v.dot(bonds[n].v);
                    R1 << cos(f*alpha), -sin(f*alpha),
                          sin(f*alpha), cos(f*alpha);
                    R2 = R1.transpose();

                    j = neighbors[i][n];

                    S1 = get_spins(i);
                    S2 = get_spins(j);

                    E0 += bonds[n].bondfunc(S1, S2);

                    E1 += bonds[n].bondfunc(S1, R1*S2);
                    Em1 += bonds[n].bondfunc(S1, R2*S2);

                    E2 += bonds[n].bondfunc(S1, R1*R1*S2);
                    Em2 += bonds[n].bondfunc(S1, R2*R2*S2);
                }
            }

            // Compute derivates from finite difference
            double dE = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha;
            double ddE = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/(alpha*alpha);
            
            return vector<double>{dE/(2.*V), ddE/(2.*V)};
        }

        inline Vector2f get_magnetization(int s) {
            Vector2f M = Vector2f::Constant(0);
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        M += get_spins(flat_idx(n1, n2, n3, s));
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
            return M/sl;
        }

        void over_relaxation_mutation(int i) {
            Vector2f H; H << 0., 0.;
            int k1; int k2; int k3; int ks;
            int j;
            for (int n = 0; n < bonds.size(); n++) {
                j = neighbors[i][n];
                H += get_spins(j);
            }

            this->mut.i = i;
            this->mut.dS = -2*get_spins(i)+ 2.*get_spins(i).dot(H)/pow(H.norm(),2) * H;
        }

        void metropolis_mutation(int i) {
            float dp = sigma*(float(r())/float(RAND_MAX) - 0.5)*2.*PI;
            Vector2f S1 = get_spins(i);
            Vector2f S2; S << cos(dp)*S1[0]
                           - sin(dp)*S1[1],
                             cos(dp)*S1[1]
                           + sin(dp)*S1[0];

            // Store mutation for consideration
            this->mut.i = i;
            this->mut.dS = S2 - S1;
        }

        void generate_mutation() {
            mut_counter++;
            mut_counter = mut_counter % V;

            if (mut_counter == 0) {
                mut_mode++;
            }

            if (mut_mode < 10) {
                over_relaxation_mutation(mut_counter);
            } else if (mut_mode < 14) {
                metropolis_mutation(mut_counter);
            } else {
                metropolis_mutation(mut_counter);
                mut_mode = 0;
            }
        }

        void accept_mutation() {
            return;
        }

        void reject_mutation() {
            set_spins(get_spins(mut.i) - mut.dS, mut.i);
        }

        virtual const float onsite_energy(int i)=0;

        virtual const float bond_energy(int i) {
            float E = 0.;
            int j;
            for (int n = 0; n < bonds.size(); n++) {
                j = neighbors[i][n];
                E += 0.5*bonds[n].bondfunc(get_spins(i), get_spins(j));
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
            set_spin(mut.i, get_spin(mut.i) + mut.dS);
            float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

            return E2 - E1;
        }

        void save_spins(string filename) {
            ofstream output_file;
            output_file.open(filename);
            output_file << N1 << "\t" << N2 << "\t" << N3 << "\t" << sl << endl;

            Vector2f S;
            Vector4i idxs;
            int n1; int n2; int n3; int s;
            for (int i = 0; i < V; i++) {
                S = get_spins(i);
                output_file << "(" << S[0] << ", " << S[1] << ")";

                n1 = idxs[0]; n2 = idxs[1]; n3 = idxs[2]; s = idxs[3];
                if (!(n1+1 == N1 && n2+1 == N2 && n3+1 == N3 && s+1 == sl)) { output_file << ","; }
            }
            output_file.close();
        }
};

template <class XYModel>
vector<float> MagnetizationLogItem(XYModel *model) {
    return vector<float>{model->get_magnetization().norm(), model->energy()};
}

#endif

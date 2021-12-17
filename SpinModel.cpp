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

struct Bond {
    int d1;
    int d2;
    int d3;
    int ds;
    Vector3f v;
    function<float(Vector2f, Vector2f)> bondfunc;
};
    
class LatticeIterator {
    public:
        int n1; int n2; int n3; int s;

        int N1; int N2; int N3; int sl;

        int counter;
        int t;

        LatticeIterator(int N1, int N2, int N3, int sl) {
            this->n1 = 0; this->n2 = 0; this->n3 = 0; this->s = 0;

            this->N1 = N1; this->N2 = N2; this->N3 = N3; this->sl = sl;

            this->counter = 0;
        }

        void next() {
            counter++;

            t = counter;
            s = t % sl;
            t = (t - s)/sl;
            n1 = t % N1;
            t = (t - n1)/N1;
            n2 = t % N2;
            t = (t - n2)/N2;
            n3 = t % N3;

            if (counter == N1*N2*N3*sl) {
                counter = 1;
            }
        }

        void reset() {
            n1 = 0; n2 = 0; n3 = 0; sl = 0;
            counter = 0;
        }
};

class SpinModel : virtual public MCModel {
    // Generic 3D Heisenberg model
    private:
        // Internal normally distributed random number generator
        GaussianDist *dist;

        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct SpinMutation {
            int n1;
            int n2;
            int n3;
            int s;
            Vector3f dS;
        };


    public:
        int sl;
        int N1;
        int N2;
        int N3;
        float acceptance;
        float sigma;
        vector<vector<vector<vector<Vector3f>>>> spins;
        vector<vector<vector<vector<vector<Vector4i>>>>> neighbors;
        vector<Bond> bonds;

        bool random_selection;
        LatticeIterator* iter;
        minstd_rand r;

        int mut_counter;
        bool mutation_mode;

        // Mutation being considered is stored as an attribute of the model
        SpinMutation mut;

        SpinModel() {}

        SpinModel(int sl, int N1, int N2 = -1, int N3 = -1) {
            this->sl = sl;
            this->N1 = N1;
            if (N2 == -1) { this->N2 = N1; } else { this->N2 = N2; }
            if (N3 == -1) { this->N3 = N1; } else { this->N3 = N3; }
            this->spins = vector<vector<vector<vector<Vector3f>>>>(this->N1,
                                 vector<vector<vector<Vector3f>>>(this->N2,
                                        vector<vector<Vector3f>>(this->N3,
                                               vector<Vector3f>(sl)))); 

            this->neighbors = vector<vector<vector<vector<vector<Vector4i>>>>>(this->N1,
                                     vector<vector<vector<vector<Vector4i>>>>(this->N2,
                                            vector<vector<vector<Vector4i>>>(this->N3,
                                                   vector<vector<Vector4i>>(this->sl,
                                                          vector<Vector4i>(0))))); 

            this->randomize_spins();

            this->acceptance = 0.5;
            this->sigma = 0.25;

            this->dist = new GaussianDist(0., 1.0);
            this->iter = new LatticeIterator(N1, N2, N3, sl);
            this->r.seed(rand());

            this->mut_counter = 0;
            this->mutation_mode = true;
            this->random_selection = false;

        }

        void randomize_spins() {
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            // For each site, initialize spin randomly
                            spins[n1][n2][n3][s] = Vector3f::Random(3);
                            // Then normalize to |S| = 1
                            spins[n1][n2][n3][s] = spins[n1][n2][n3][s].normalized();
                        }
                    }
                }
            }
        }

        void add_bond(Bond b) {
            Vector4i v;
            this->bonds.push_back(b);
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            v << mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3), mod(s + b.ds, sl);
                            neighbors[n1][n2][n3][s].push_back(v);
                        }
                    }
                }
            }
        }

        inline vector<double> twist_stiffness() {
            // Returns the first and second derivative in response to a phase twist
            float alpha = 0.01;
            float f;
            Matrix3f R1; 
            Matrix3f R2;

            double E0 = 0.;
            double E1 = 0.;
            double E2 = 0.;
            double Em1 = 0.;
            double Em2 = 0.;


            Vector3f S1;
            Vector3f S2;
            int k1; int k2; int k3; int ks;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            for (int n = 0; n < bonds.size(); n++) {
                                f = bonds[0].v.dot(bonds[n].v);
                                R1 << cos(f*alpha), -sin(f*alpha), 0,
                                      sin(f*alpha), cos(f*alpha), 0.,
                                      0., 0., 1.;
                                R2 = R1.transpose();

                                k1 = neighbors[n1][n2][n3][s][n][0]; k2 = neighbors[n1][n2][n3][s][n][1]; 
                                k3 = neighbors[n1][n2][n3][s][n][2]; ks = neighbors[n1][n2][n3][s][n][3];

                                S1 = spins[n1][n2][n3][s];
                                S2 = spins[k1][k2][k3][ks];

                                E0 += bonds[n].bondfunc(S1, S2);

                                E1 += bonds[n].bondfunc(S1, R1*S2);
                                Em1 += bonds[n].bondfunc(S1, R2*S2);

                                E2 += bonds[n].bondfunc(S1, R1*R1*S2);
                                Em2 += bonds[n].bondfunc(S1, R2*R2*S2);
                            }
                        }
                    }
                }
            }

            // Compute derivates from finite difference
            double dE = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha;
            double ddE = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/(alpha*alpha);
            
            //return vector<double>{dE, ddE};
            return vector<double>{dE/2., ddE/2.};
        }

        inline Vector3f get_magnetization(int s) {
            Vector3f M = Vector3f::Constant(0);
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        M += spins[n1][n2][n3][s];
                    }
                }
            }
            
            return M/(N1*N2*N3);
        }

        inline Vector3f get_magnetization() {
            Vector3f M = Vector3f::Constant(0);
            for (int s = 0; s < sl; s++) {
                M += get_magnetization(s);
            }
            return M/(sl);
        }

        void over_relaxation_mutation(int n1, int n2, int n3, int s) {
            Vector3f H; H << 0., 0., 0.;
            int k1; int k2; int k3; int ks;
            for (int n = 0; n < bonds.size(); n++) {
                k1 = neighbors[n1][n2][n3][s][n][0]; k2 = neighbors[n1][n2][n3][s][n][1]; 
                k3 = neighbors[n1][n2][n3][s][n][2]; ks = neighbors[n1][n2][n3][s][n][3];
                H += this->spins[k1][k2][k3][ks];
            }

            this->mut.n1 = n1;
            this->mut.n2 = n2;
            this->mut.n3 = n3;
            this->mut.s = s;
            this->mut.dS = -2*this->spins[n1][n2][n3][s] + 2.*this->spins[n1][n2][n3][s].dot(H)/pow(H.norm(),2) * H;
        }

        void metropolis_mutation(int n1, int n2, int n3, int s) {
            if (acceptance > 0.5) {
                sigma = min(2., 1.01*sigma);

            } else {
                sigma = max(0.05, 0.99*sigma);
            }

            // Randomly generate mutation
            Vector3f Gamma;
            Gamma << dist->sample(), dist->sample(), dist->sample();
            Vector3f S2 = (spins[n1][n2][n3][s] + this->sigma*Gamma).normalized();


            // Store mutation for consideration
            this->mut.n1 = n1;
            this->mut.n2 = n2;
            this->mut.n3 = n3;
            this->mut.s = s;
            this->mut.dS = S2 - spins[n1][n2][n3][s];
        }

        void generate_mutation() {
            mut_counter++;
            if ((mut_counter % (N1*N2*N3*sl))%5 == 0) { mutation_mode = false; mut_counter = 1; } else { mutation_mode = true; }

            int n1; int n2; int n3; int s;
            if (random_selection) {
                n1 = r() % N1;
                n2 = r() % N2;
                if (this->N3 == 1) { n3 = 0; } else { n3 = r() % N3; }

                if (sl == 1) { s = 0; } else { s = r() % sl; }
            } else {
                n1 = iter->n1;
                n2 = iter->n2;
                n3 = iter->n3;
                s = iter->s;
                iter->next();
            }


            if (mutation_mode) {
                over_relaxation_mutation(n1, n2, n3, s);
            } else {
                metropolis_mutation(n1, n2, n3, s);
            }
        }

        void accept_mutation() {
            return;
        }

        void reject_mutation() {
            this->spins[mut.n1][mut.n2][mut.n3][mut.s] -= mut.dS;
        }

        virtual const float onsite_energy(int n1, int n2, int n3, int s)=0;

        virtual const float bond_energy(int n1, int n2, int n3, int s) {
            float E = 0.;
            int k1; int k2; int k3; int ks;
            for (int n = 0; n < bonds.size(); n++) {
                k1 = neighbors[n1][n2][n3][s][n][0]; k2 = neighbors[n1][n2][n3][s][n][1]; 
                k3 = neighbors[n1][n2][n3][s][n][2]; ks = neighbors[n1][n2][n3][s][n][3];
                E += 0.5*bonds[n].bondfunc(spins[n1][n2][n3][s], spins[k1][k2][k3][ks]);
            }

            return E;
        }

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
            this->spins[mut.n1][mut.n2][mut.n3][mut.s] += mut.dS;
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
                            output_file << "[" << spins[n1][n2][n3][s][0] << " " 
                                               << spins[n1][n2][n3][s][1] << " " 
                                               << spins[n1][n2][n3][s][2] << "]";
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

template <class SpinModel>
vector<float> MagnetizationLogItem(SpinModel *model) {
    return vector<float>{model->get_magnetization().norm(), model->energy()};
}

#endif

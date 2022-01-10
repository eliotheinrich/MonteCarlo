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
    function<float(Vector3f, Vector3f)> bondfunc;
};    

class SpinModel : virtual public MCModel {
    // Generic 3D Heisenberg model
    private:
        // Internal normally distributed random number generator
        GaussianDist *dist;

        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct SpinMutation {
            int i;
            Vector3f dS;
        };


    public:
        int sl;
        int N1;
        int N2;
        int N3;
        int V;

        float acceptance;
        float sigma;
        vector<Vector3f> spins;
        vector<vector<int>> neighbors;
        vector<Bond> bonds;

        minstd_rand r;

        int mut_counter;
        int mut_mode;

        // Mutation being considered is stored as an attribute of the model
        SpinMutation mut;

        SpinModel() {}

        SpinModel(int sl, int N1, int N2 = -1, int N3 = -1) {
            this->sl = sl;
            this->N1 = N1;
            if (N2 == -1) { this->N2 = N1; } else { this->N2 = N2; }
            if (N3 == -1) { this->N3 = N1; } else { this->N3 = N3; }
            this->V = N1*N2*N3*sl;

            this->spins = vector<Vector3f>(V);
            this->neighbors = vector<vector<int>>(V, vector<int>(0));

            this->randomize_spins();

            this->acceptance = 0.5;
            this->sigma = 0.25;

            this->dist = new GaussianDist(0., 1.0);
            this->r.seed(rand());

            this->mut_counter = 0;
        }

        inline int flat_idx(int n1, int n2, int n3, int s) {
            return n1 + N1*(n2 + N2*(n3 + N3*s));
        }

        inline Vector4i tensor_idx(int i) {
            int n1 = i % N1;
            i = i / N1;
            int n2 = i % N2;
            i = i / N2;
            int n3 = i % N3;
            i = i / N3;
            int s = i % sl;
            Vector4i v; v << n1, n2, n3, s;
            return v;
        }

        inline void set_spin(int n1, int n2, int n3, int s, Vector3f v) {
            spins[flat_idx(n1, n2, n3, s)] = v;
        }

        inline void set_spin(int i, Vector3f v) {
            spins[i] = v;
        }

        inline const Vector3f get_spin(int n1, int n2, int n3, int s) {
            return spins[flat_idx(n1, n2, n3, s)];
        }

        inline const Vector3f get_spin(int i) {
            return spins[i];
        }

        void randomize_spins() {
            for (int i = 0; i < V; i++) {
                set_spin(i, Vector3f::Random(3).normalized());
            }
        }

        void add_bond(Bond b) {
            this->bonds.push_back(b);
            int idx; int neighbor_idx;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            idx = flat_idx(n1, n2, n3, s);
                            neighbor_idx = flat_idx(mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3), mod(s + b.ds, sl));
                            neighbors[idx].push_back(neighbor_idx);
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
            int j;
            for (int i = 0; i < V; i++) {
                for (int n = 0; n < bonds.size(); n++) {
                    f = bonds[0].v.dot(bonds[n].v);
                    R1 << cos(f*alpha), -sin(f*alpha), 0,
                          sin(f*alpha), cos(f*alpha), 0.,
                          0., 0., 1.;
                    R2 = R1.transpose();

                    j = neighbors[i][n];

                    S1 = get_spin(i);
                    S2 = get_spin(j);

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
            
            return vector<double>{dE/(2.*sqrt(V)), ddE/(2*V)};
        }

        inline Vector3f get_magnetization() {
            Vector3f M = Vector3f::Constant(0);
            for (int i = 0; i < V; i++) {
                M += spins[i];
            }
            
            return M/(N1*N2*N3*sl);
        }

        vector<float> correlation_function(int i, int a = 2, int b = 2) {
            vector<float> Cij = vector<float>(V); 

            int j;
            Vector4i idxs = tensor_idx(i);
            int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            j = flat_idx(n1, n2, n3, s);
                            Cij[j] = spins[flat_idx((m1 + n1)%N1, 
                                                    (m2 + n2)%N2, 
                                                    (m3 + n3)%N3, 
                                                    (s + k)%sl)][a]*spins[i][b];
                        }
                    }
                }
            }
            return Cij;
        }

        vector<float> full_correlation_function(int i) {
            vector<float> Cij = vector<float>(V); 

            int j;
            Vector4i idxs = tensor_idx(i);
            int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            j = flat_idx(n1, n2, n3, s);
                            Cij[j] = spins[flat_idx((m1 + n1)%N1, 
                                                    (m2 + n2)%N2, 
                                                    (m3 + n3)%N3, 
                                                    (s + k)%sl)].dot(spins[i]);
                        }
                    }
                }
            }
            return Cij;
        }

        float skyrmion_density(int i) {
            int j;

            Vector3f dSdX; dSdX << 0., 0., 0.;
            Vector3f dSdY; dSdY << 0., 0., 0.;
            for (int n = 0; n < bonds.size(); n++) {
                j = neighbors[i][n];
                if (bonds[n].v[0] != 0.) {
                    dSdX += bonds[n].v[0]*(spins[j] - spins[i]);
                }
                if (bonds[n].v[1] != 0.) {
                    dSdX += bonds[n].v[1]*(spins[j] - spins[i]);
                }

            }
            dSdX = dSdX/bonds.size();
            dSdY = dSdY/bonds.size();

            return spins[i].dot(dSdX.cross(dSdY));
        }

        vector<float> skyrmion_correlation_function(int i) {
            vector<float> Cij = vector<float>(V); 

            int j;
            Vector4i idxs = tensor_idx(i);
            int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            j = flat_idx(n1, n2, n3, s);
                            Cij[j] = skyrmion_density(flat_idx((m1 + n1)%N1, 
                                                               (m2 + n2)%N2, 
                                                               (m3 + n3)%N3, 
                                                               (s + k)%sl))*skyrmion_density(i);
                        }
                    }
                }
            }
            return Cij;
        }

        void over_relaxation_mutation(int i) {
            Vector3f H; H << 0., 0., 0.;
            int j;
            for (int n = 0; n < bonds.size(); n++) {
                j = neighbors[i][n];
                H += spins[j];
            }

            this->mut.i = i;
            this->mut.dS = -2*spins[i] + 2.*spins[i].dot(H)/pow(H.norm(),2) * H;
        }

        void metropolis_mutation(int i) {
            if (acceptance > 0.5) {
                sigma = min(2., 1.01*sigma);

            } else {
                sigma = max(0.05, 0.99*sigma);
            }

            // Randomly generate mutation
            Vector3f Gamma;
            Gamma << dist->sample(), dist->sample(), dist->sample();
            Vector3f S2 = (spins[i] + this->sigma*Gamma).normalized();


            // Store mutation for consideration
            this->mut.i = i;
            this->mut.dS = S2 - spins[i];
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
            spins[mut.i] = spins[mut.i] - mut.dS;
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
            spins[mut.i] = spins[mut.i] + mut.dS;
            float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

            return E2 - E1;
        }

        // Saves current spin configuration
        void save_spins(string filename) {
            ofstream output_file;
            output_file.open(filename);
            output_file << sl << "\t" << N1 << "\t" << N2 << "\t" << N3 << endl;
            int i; Vector3f S;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            i = flat_idx(n1, n2, n3, s);
                            S = get_spin(i);
                            output_file << "[" << S[0] << " " 
                                               << S[1] << " " 
                                               << S[2] << "]";
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

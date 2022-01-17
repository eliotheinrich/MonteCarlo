#ifndef SPINMC_
#define SPINMC_

#include <iostream>
#include <functional>
#include <vector>
#include <unordered_set>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <Eigen/Dense>
#include <fstream>
#include "MonteCarlo.cpp"
#include "Utility.cpp"

    
class XYModel : virtual public MCModel {
    // Generic 3D XY model
    private:
        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct XYMutation {
            int i;
            Eigen::Vector2f dS;
        };

        struct XYBond {
            int d1;
            int d2;
            int d3;
            int ds;
            Eigen::Vector3f v;
            std::function<float(Eigen::Vector2f, Eigen::Vector2f)> bondfunc;
        };


    public:
        int sl;
        int N1;
        int N2;
        int N3;
        int V;

        float acceptance;
        float sigma;
        std::vector<Eigen::Vector2f> spins;
        std::vector<std::vector<int>> neighbors;
        std::vector<XYBond> bonds;

        std::unordered_set<int> s;

        std::minstd_rand r;

        bool cluster;
        float dE;

        // Mutation being considered is stored as an attribute of the model
        XYMutation mut;

        XYModel() {}

        XYModel(int sl, int N1, int N2 = -1, int N3 = -1) {
            this->sl = sl;
            this->N1 = N1;
            if (N2 == -1) { this->N2 = N1; } else { this->N2 = N2; }
            if (N3 == -1) { this->N3 = N1; } else { this->N3 = N3; }
            this->V = N1*N2*N3*sl;

            this->spins = std::vector<Eigen::Vector2f>(V);
            this->neighbors = std::vector<std::vector<int>>(V, std::vector<int>(0));

            this->randomize_spins();

            this->acceptance = 0.5;
            this->sigma = 0.25;

            this->r.seed(rand());

            this->mut.i = 0;
            this->cluster = true;
        }

        inline const int flat_idx(int n1, int n2, int n3, int s) {
            return n1 + N1*(n2 + N2*(n3 + N3*s));
        }

        inline const Eigen::Vector4i tensor_idx(int i) {
            int n1 = i % N1;
            i = i / N1;
            int n2 = i % N2;
            i = i / N2;
            int n3 = i % N3;
            i = i / N3;
            int s = i % sl;
            Eigen::Vector4i v; v << n1, n2, n3, s;
            return v;
        }

        inline void set_spin(int n1, int n2, int n3, int s, Eigen::Vector2f v) {
            spins[flat_idx(n1, n2, n3, s)] = v;
        }


        inline const Eigen::Vector2f get_spin(int n1, int n2, int n3, int s) {
            return spins[flat_idx(n1, n2, n3, s)];
        }

        void randomize_spins() {
            float p;
            Eigen::Vector2f v;

            for (int i = 0; i < V; i++) {
                p = 2*PI*float(r())/float(RAND_MAX);
                spins[i] << cos(p), sin(p);
            }
        }

        void add_bond(int d1, int d2, int d3, int ds, Eigen::Vector3f v, std::function<float(Eigen::Vector2f, Eigen::Vector2f)> bondfunc) {
            XYBond b{d1, d2, d3, ds, v, bondfunc};
            this->bonds.push_back(b);
            int i; int j;
            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        for (int s = 0; s < sl; s++) {
                            i = flat_idx(n1, n2, n3, s);
                            j = flat_idx(mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3), mod(s + b.ds, sl));
                            neighbors[i].push_back(j);
                        }
                    }
                }
            }
        }

        inline std::vector<double> twist_stiffness() {
            // Returns the first and second derivative in response to a phase twist
            float alpha = 0.01;
            float f;
            Eigen::Matrix2f R1; 
            Eigen::Matrix2f R2;
            Eigen::Matrix2f R3;

            double E0 = 0.;
            double E1 = 0.;
            double E2 = 0.;
            double E3 = 0.;
            double Em1 = 0.;
            double Em2 = 0.;
            double Em3 = 0.;


            Eigen::Vector2f S1;
            Eigen::Vector2f S2;
            int j;
            for (int i = 0; i < V; i++) {
                for (int n = 0; n < bonds.size(); n++) {
                    f = bonds[0].v.dot(bonds[n].v);
                    R1 << std::cos(f*alpha), -std::sin(f*alpha),
                          std::sin(f*alpha), std::cos(f*alpha);
                    R2 = R1*R1;
                    R3 = R1*R1*R1;

                    j = neighbors[i][n];

                    S1 = spins[i];
                    S2 = spins[j];

                    E0 += bonds[n].bondfunc(S1, S2);

                    E1 += bonds[n].bondfunc(S1, R1*S2);
                    Em1 += bonds[n].bondfunc(S1, R1.transpose()*S2);

                    E2 += bonds[n].bondfunc(S1, R2*S2);
                    Em2 += bonds[n].bondfunc(S1, R2.transpose()*S2);

                    E3 += bonds[n].bondfunc(S1, R3*S2);
                    Em3 += bonds[n].bondfunc(S1, R3.transpose()*S2);
                }
            }

            // Compute derivates from finite difference
            double d1E = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha;
            double d2E = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/pow(alpha, 2);
            double d3E = (1./8.*Em3 - 1.*Em2 + 13./8.*Em1 - 13./8.*E1 + 1.*E2 - 1./8.*E3)/pow(alpha, 3);
            double d4E = (-1./6.*Em3 + 2.*Em2 - 13./2.*Em1 + 28./3.*E0 - 13./2.*E1 + 2.*E2 - 1./6.*E3)/pow(alpha, 4);
            
            return std::vector<double>{d1E/2., d2E/2., d3E/2., d4E/2.};
        }

        inline Eigen::Vector2f get_magnetization() {
            Eigen::Vector2f M = Eigen::Vector2f::Constant(0);
            for (int i = 0; i < V; i++) {
                M += spins[i];
            }
            
            return M/(N1*N2*N3*sl);
        }

        void cluster_update() {
            s.clear();

            float p = (float) r()/RAND_MAX;
            Eigen::Vector2f ax; ax << std::cos(p), std::sin(p);
            ax = ax.normalized();

            int i = r() % V;
            spins[i] = spins[i] - 2*spins[i].dot(ax)*ax;

            std::stack<int> c;
            c.push(i);

            float dE;
            Eigen::Vector2f new_S;
            int m; int j;
            while (!c.empty()) {
                m = c.top();
                c.pop();

                // Mark m as visited
                s.insert(m);

                for (int n = 0; n < bonds.size(); n++) {
                    j = neighbors[m][n];
                    // Check if each neighbor has been visited
                    if (!s.count(j)) {
                        // With appropriate probability, add neighbor to stack and flip it
                        new_S = spins[j] - 2*spins[j].dot(ax)*ax;
                        dE = bonds[n].bondfunc(spins[m], new_S) - bonds[n].bondfunc(spins[m], spins[j]);
                        if ((float) r()/RAND_MAX < 1. - exp(dE/T)) {
                            c.push(j);
                            spins[j] = new_S;
                        }
                    }
                }
            }
        }

        void metropolis_mutation() {
            float dp = sigma*(float(r())/float(RAND_MAX) - 0.5)*2.*PI;
            Eigen::Vector2f S1 = spins[mut.i];
            Eigen::Vector2f S2; S2 << std::cos(dp)*S1[0]
                             - std::sin(dp)*S1[1],
                               std::cos(dp)*S1[1]
                             + std::sin(dp)*S1[0];

            // Store mutation for consideration
            this->mut.dS = S2 - S1;
        }

        void generate_mutation() {
            if (cluster) {
                cluster_update();
            } else {
                mut.i++;
                if (mut.i == V) {
                    mut.i = 0;
                }
                metropolis_mutation();
            }
        }

        void accept_mutation() {
            return;
        }

        void reject_mutation() {
            if (!cluster) {
                spins[mut.i] = spins[mut.i] - mut.dS;
            }
        }

        virtual const float onsite_func(const Eigen::Vector2f& S)=0;

        virtual const float onsite_energy(int i) {
            return onsite_func(spins[i]);
        }

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
            if (cluster) {
                return -1.;
            } else {
                float E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
                spins[mut.i] = spins[mut.i] + mut.dS;
                float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

                return E2 - E1;
            }
        }

        void save_spins(std::string filename) {
            std::ofstream output_file;
            output_file.open(filename);
            output_file << N1 << "\t" << N2 << "\t" << N3 << "\t" << sl << "\n";

            Eigen::Vector2f S;
            for (int i = 0; i < V; i++) {
                S = spins[i];
                output_file << S[0] << "\t" << S[1];
                if (i < V-1) { output_file << "\t"; }
            }
            output_file.close();
        }
};

#endif

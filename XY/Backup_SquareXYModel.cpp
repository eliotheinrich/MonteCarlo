#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../XYModel.cpp"
#include "../Utility.cpp"


using namespace std;
using namespace Eigen;

class SquareXYModel : public XYModel {
    public:
        int N;
        int L;
        float J;
        float B;
        float Bp;
        float Bx;
        float By;

        int counter;
        bool mutation_mode;


        SquareXYModel(int N, int L, float J, float B, float Bp) : XYModel(1, N, N, L) {
            this->N = N;
            this->L = L;
            this->J = J;
            this->B = B;
            this->Bp = Bp;
            this->Bx = B*cos(Bp);
            this->By = B*sin(Bp);

            function<float(Vector2f, Vector2f)> dotfunc = [J](Vector2f S1, Vector2f S2) {
                return -J*S1.dot(S2);
            };

            this->add_bond(Bond{1,0,0,0,dotfunc});
            this->add_bond(Bond{-1,0,0,0,dotfunc});
            this->add_bond(Bond{0,1,0,0,dotfunc});
            this->add_bond(Bond{0,-1,0,0,dotfunc});

            this->counter = 0;
            this->mutation_mode = true;
        }

        const float onsite_energy(int n1, int n2, int n3, int s) {
            float E = 0;

            // Onsite interactions
            E -= this->Bx*this->spins[n1][n2][n3][s][0] + this->By*this->spins[n1][n2][n3][s][1];

            return E;
        }

        const float bond_energy(int n1, int n2, int n3, int s) {
            float E = 0;

            // NN interactions

            // PBC boundary conditions
            E -= this->spins[n1][n2][n3][s].dot(this->spins[(n1+1)%N][n2][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[mod(n1-1,N)][n2][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[n1][(n2+1)%N][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[n1][mod(n2-1,N)][n3][s]);

            return 0.5*J*E;
        }
/*
        inline float spin_stiffness(float T) {
            float alpha = 0.01;
            Matrix2f R1; R1 << cos(alpha), -sin(alpha),
                               sin(alpha), cos(alpha);
            Matrix2f R2 = R1.transpose();

            double E0 = 0.;
            double E1 = 0.;
            double E2 = 0.;
            double Em1 = 0.;
            double Em2 = 0.;


            for (int n1 = 0; n1 < N1; n1++) {
                for (int n2 = 0; n2 < N2; n2++) {
                    for (int n3 = 0; n3 < N3; n3++) {
                        E0 += -J*spins[n1][n2][n3][0].dot(spins[(n1+1)%N1][n2][n3][0]);

                        E1 += -J*spins[n1][n2][n3][0].dot(R1*spins[(n1+1)%N1][n2][n3][0]);
                        Em1 += -J*spins[n1][n2][n3][0].dot(R2*spins[(n1+1)%N1][n2][n3][0]);

                        E2 += -J*spins[n1][n2][n3][0].dot(R1*R1*spins[(n1+1)%N1][n2][n3][0]);
                        Em2 += -J*spins[n1][n2][n3][0].dot(R2*R2*spins[(n1+1)%N1][n2][n3][0]);
                    }
                }
            }

            double dE = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha;
            double ddE = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/(alpha*alpha);

            return (ddE - (1./T)*dE*dE)/(N*N);
        }

        void over_relaxation_mutation() {
            int n1 = XYModel::iter->n1;
            int n2 = XYModel::iter->n2;
            int n3 = XYModel::iter->n3;
            int s = XYModel::iter->s;
            XYModel::iter->next();

            Vector2f H; H << 0., 0.;
            H += this->spins[(n1+1)%N][n2][n3][0];
            H += this->spins[mod(n1-1,N)][n2][n3][0];
            H += this->spins[n1][(n2+1)%N][n3][0];
            H += this->spins[n1][mod(n2-1,N)][n3][0];

            this->mut.n1 = n1;
            this->mut.n2 = n2;
            this->mut.n3 = n3;
            this->mut.dS = -2*this->spins[n1][n2][n3][0] + 2.*this->spins[n1][n2][n3][0].dot(H)/pow(H.norm(),2) * H;
        }

        void generate_mutation() {
            counter++;
            if ((counter % (N*N*L))%5 == 0) { mutation_mode = false; } else { mutation_mode = true; }

            if (mutation_mode) {
                over_relaxation_mutation();
            } else {
                XYModel::generate_mutation();
            }
        }
        */

};


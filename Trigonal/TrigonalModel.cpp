#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../SpinModel.cpp"
#include "../Utility.cpp"


using namespace std;
using namespace Eigen;


class TrigonalModel : public SpinModel {
    public:
        int N;
        int L;
        float J1;
        float J2;
        float K1;
        float K2;
        float K3;

        Vector3f B; 
        float S;

        Matrix3f R;
        int counter;


        TrigonalModel(int N, int L, float J1, float J2, float K1, float K2, float K3,
                                    Vector3f B, float S) : SpinModel(1, N, N, L) {

            this->N = N;
            this->L = L;
            this->J1 = J1;
            this->J2 = J2;
            this->K1 = K1;
            this->K2 = K2;
            this->K3 = K3;
            this->B = B;
            this->S = S;
            cout << "J1 = " << this->J1 << endl;
            cout << "J2 = " << this->J2 << endl;
            cout << "K1 = " << this->K1 << endl;
            cout << "K2 = " << this->K2 << endl;
            cout << "K3 = " << this->K3 << endl;
            cout << "S = " << this->S << endl;

            this->R << sqrt(3)/2., -0.5, 0.,
                       0.5, sqrt(3)/2., 0.,
                       0., 0., 1.;
            this->counter = 0;

            function<float(Vector3f S1, Vector3f S2)> bondfunc = [J1](Vector3f S1, Vector3f S2) {
                return -J1*S1.dot(S2);
            };

            this->add_bond(Bond{1,0,0,0, bondfunc});
            this->add_bond(Bond{-1,0,0,0, bondfunc});
            this->add_bond(Bond{0,1,0,0, bondfunc});
            this->add_bond(Bond{0,-1,0,0, bondfunc});
            this->add_bond(Bond{1,-1,0,0, bondfunc});
            this->add_bond(Bond{-1,1,0,0, bondfunc});
        }

/*
        void generate_mutation() {
            counter++;
            if (counter % 10 == 0) {
                // Randomly select a site to mutate
                int rand_n1 = rand() % N1;
                int rand_n2 = rand() % N2;
                int rand_n3;
                if (this->N3 == 1) { rand_n3 = 0; } else { rand_n3 = rand() % N3; }

                int rand_s;
                if (sl == 1) { rand_s = 0; } else { rand_s = rand() % sl; }

                // Randomly generate mutation
                Vector3f S2;
                if (rand() % 2 == 0) {
                    S2 = R*spins[rand_n1][rand_n2][rand_n3][rand_s];
                } else {
                    S2 = (R.transpose())*spins[rand_n1][rand_n2][rand_n3][rand_s];
                }
                

                // Store mutation for consideration
                this->mut.n1 = rand_n1;
                this->mut.n2 = rand_n2;
                this->mut.n3 = rand_n3;
                this->mut.s = rand_s;
                this->mut.dS = S2 - spins[rand_n1][rand_n2][rand_n3][rand_s];
            } else {
                SpinModel::generate_mutation();
            }   
        }
*/
        const float onsite_energy(int n1, int n2, int n3, int s) {
            float E = 0;

            // Onsite interactions
            E -= S*B.dot(this->spins[n1][n2][n3][s]);
            float x = this->spins[n1][n2][n3][s][0];
            float y = this->spins[n1][n2][n3][s][1];
            float z = this->spins[n1][n2][n3][s][2];
            E += K1*S*S*z*z;
            E += K2*S*pow(x*x+y*y,2);
            E += K3*S*cos(6*atan2(y, x))*pow(x*x+y*y,3); // Sixfold magnetocrystalline field

            return E;
        }

        const float bond_energy(int n1, int n2, int n3, int s) {
            float E = 0;

            // NN interactions
            E -= this->spins[n1][n2][n3][s].dot(this->spins[n1][(n2+1)%N][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[n1][mod(n2-1, N)][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[(n1+1)%N][n2][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[mod(n1-1, N)][n2][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[(n1+1)%N][mod(n2-1, N)][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[mod(n1-1, N)][(n2+1)%N][n3][s]);
            E *= 0.5*J1;

            // PBC on interlayer coupling
            //E += 0.5*J2*S*S*this->spins[n1][n2][n3][s].dot(this->spins[n1][n2][mod(n3-1, L)][s]);
            //E += 0.5*J2*S*S*this->spins[n1][n2][n3][s].dot(this->spins[n1][n2][(n3+1)% L][s]);

            // OBC on interlayer coupling 
            //if (s > 0) {
            //    E += 0.5*J2*S*S*this->spins[n1][n2][n3][s].dot(this->spins[n1][n2][n3-1][s]);
            //} 
            //if (s < L-1) {
            //    E += 0.5*J2*S*S*this->spins[n1][n2][n3][s].dot(this->spins[n1][n2][n3+1][s]);
            //}

            return E;
        }
};


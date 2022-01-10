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

        Matrix3f R;
        int mut_type;


        TrigonalModel(int N, int L, float J1, float J2, float K1, float K2, float K3,
                                    Vector3f B) : SpinModel(1, N, N, L) {

            this->N = N;
            this->L = L;
            this->J1 = J1;
            this->J2 = J2;
            this->K1 = K1;
            this->K2 = K2;
            this->K3 = K3;
            this->B = B;


            if (false) {
                cout << "J1 = " << J1 << endl;
                cout << "J2 = " << J2 << endl;
                cout << "K1 = " << K1 << endl;
                cout << "K2 = " << K2 << endl;
                cout << "K3 = " << K3 << endl;

                cout << "Bx = " << B[0] << endl;
                cout << "By = " << B[1] << endl;
                cout << "Bz = " << B[2] << endl;
            }

            this->R << sqrt(3)/2., -0.5, 0.,
                       0.5, sqrt(3)/2., 0.,
                       0., 0., 1.;
            this->mut_type = 0;

            function<float(Vector3f S1, Vector3f S2)> bondfunc = [J1](Vector3f S1, Vector3f S2) {
                return -J1*S1.dot(S2);
            };

            Vector3f v1; v1 << 1., 0., 0.;
            Vector3f v2; v2 << 0.5, sqrt(3)/2., 0.;
            Vector3f v3; v3 << 0.5, -sqrt(3)/2., 0.;
            this->add_bond(Bond{1,0,0,0, v1, bondfunc});
            this->add_bond(Bond{-1,0,0,0, -v1, bondfunc});
            this->add_bond(Bond{0,1,0,0, v2, bondfunc});
            this->add_bond(Bond{0,-1,0,0, -v2, bondfunc});
            this->add_bond(Bond{1,-1,0,0, v3, bondfunc});
            this->add_bond(Bond{-1,1,0,0, -v3, bondfunc});

            function<float(Vector3f S1, Vector3f S2)> bondfunc_inter = [J2](Vector3f S1, Vector3f S2) {
                return J2*S1.dot(S2);
            };

            Vector3f v4; v4 << 0., 0., 1.;
            this->add_bond(Bond{0,0,1,0, v4, bondfunc_inter});
            this->add_bond(Bond{0,0,-1,0, -v4, bondfunc_inter});
        }

        TrigonalModel* clone() {
            TrigonalModel* new_model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
            for (int i = 0; i < V; i++) {
                new_model->set_spin(i, this->get_spin(i));
            }
            return new_model;
        }

        void over_relaxation_mutation(int i) {
            Vector3f H; H << 0., 0., 0.;
            int j;
            for (int n = 0; n < 6; n++) {
                j = neighbors[i][n];
                H -= J1*get_spin(j);
            }
            for (int n = 6; n < 8; n++) {
                j = neighbors[i][n];
                H -= J1*get_spin(j);
            }

            this->mut.i = i;
            this->mut.dS = -2*get_spin(i) + 2.*get_spin(i).dot(H)/pow(H.norm(),2) * H;
        }

        void rotation_mutation(int i) {
            Vector3f S2;
            if (r() % 2) {
                S2 = R*get_spin(i);
            } else {
                S2 = R.transpose()*get_spin(i);
            }

            this->mut.i = i;
            this->mut.dS = S2 - get_spin(i);
        }

        const float onsite_energy(int i) {
            float E = 0;

            // Onsite interactions
            Vector3f S = get_spin(i);
            E -= B.dot(S);

            E += K1*S[2]*S[2];
            //E += K2*pow(S[0]*S[0]+S[1]*S[1],2);
            E += K3*cos(6*atan2(S[1], S[0]))*pow(S[0]*S[0]+S[1]*S[1],3); // Sixfold magnetocrystalline field

            return E;
        }

        const float bond_energy(int i) {
            float E = 0;

            Vector4i idxs = tensor_idx(i);
            Vector3f S = get_spin(i);

            int n1 = idxs[0]; int n2 = idxs[1]; int n3 = idxs[2]; int s = idxs[3];

            // NN interactions
            E -= S.dot(get_spin(n1, mod(n2+1, N), n3, s));
            E -= S.dot(get_spin(n1, mod(n2-1, N), n3, s));
            E -= S.dot(get_spin(mod(n1+1, N), n2, n3, s));
            E -= S.dot(get_spin(mod(n1-1, N), n2, n3, s));
            E -= S.dot(get_spin(mod(n1+1, N), mod(n2-1, N), n3, s));
            E -= S.dot(get_spin(mod(n1-1, N), mod(n2+1, N), n3, s));
            E = 0.5*J1*E;

            // PBC on interlayer coupling
            E += 0.5*J2*S.dot(get_spin(n1, n2, mod(n3+1, L), s));
            E += 0.5*J2*S.dot(get_spin(n1, n2, mod(n3-1, L), s));

            return E;
        }

};


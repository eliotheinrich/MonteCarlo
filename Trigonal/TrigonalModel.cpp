#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../SpinModel.cpp"
#include "../Utility.cpp"

class TrigonalModel : public SpinModel {
    public:
        int N;
        int L;
        float J1;
        float J2;
        float K1;
        float K2;
        float K3;
        Eigen::Vector3f B; 

        Eigen::Matrix3f R;
        int mut_counter;
        int mut_type;


        TrigonalModel(int N, int L, float J1, float J2, float K1, float K2, float K3,
                                    Eigen::Vector3f B) : SpinModel(1, N, N, L) {

            this->N = N;
            this->L = L;
            this->J1 = J1;
            this->J2 = J2;
            this->K1 = K1;
            this->K2 = K2;
            this->K3 = K3;
            this->B = B;


            if (false) {
                std::cout << "J1 = " << J1 << std::endl;
                std::cout << "J2 = " << J2 << std::endl;
                std::cout << "K1 = " << K1 << std::endl;
                std::cout << "K2 = " << K2 << std::endl;
                std::cout << "K3 = " << K3 << std::endl;

                std::cout << "Bx = " << B[0] << std::endl;
                std::cout << "By = " << B[1] << std::endl;
                std::cout << "Bz = " << B[2] << std::endl;
            }

            this->R << sqrt(3)/2., -0.5, 0.,
                       0.5, sqrt(3)/2., 0.,
                       0., 0., 1.;
            this->mut_type = 0;
            this->mut_counter = 0; 

            std::function<float(const Eigen::Vector3f &, const Eigen::Vector3f &)> bondfunc = 
            [J1](const Eigen::Vector3f &S1, const Eigen::Vector3f &S2) {
                return -J1*S1.dot(S2);
            };

            Eigen::Vector3f v1; v1 << 1., 0., 0.;
            Eigen::Vector3f v2; v2 << 0.5, std::sqrt(3)/2., 0.;
            Eigen::Vector3f v3; v3 << 0.5, -std::sqrt(3)/2., 0.;
            this->add_bond(1,0,0,0, v1, bondfunc);
            this->add_bond(-1,0,0,0, -v1, bondfunc);
            this->add_bond(0,1,0,0, v2, bondfunc);
            this->add_bond(0,-1,0,0, -v2, bondfunc);
            this->add_bond(1,-1,0,0, v3, bondfunc);
            this->add_bond(-1,1,0,0, -v3, bondfunc);

            std::function<float(const Eigen::Vector3f &, const Eigen::Vector3f &)> bondfunc_inter = 
            [J2](const Eigen::Vector3f &S1, const Eigen::Vector3f &S2) {
                return J2*S1.dot(S2);
            };

            if ((std::abs(J2) > 1e-6) && (L > 1)) {
                Eigen::Vector3f v4; v4 << 0., 0., 1.;
                this->add_bond(0,0,1,0, v4, bondfunc_inter);
                this->add_bond(0,0,-1,0, -v4, bondfunc_inter);
            }
        }

        TrigonalModel* clone() {
            TrigonalModel* new_model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
            for (int i = 0; i < V; i++) {
                new_model->spins[i] = this->spins[i];
            }
            return new_model;
        }

        void over_relaxation_mutation() {
            Eigen::Vector3f H; H << 0., 0., 0.;
            H = H + B;

            int j;
            for (int n = 0; n < 6; n++) {
                j = neighbors[mut.i][n];
                H -= J1*spins[j];
            }
            if (bonds.size() > 7) {
                for (int n = 6; n < 8; n++) {
                    j = neighbors[mut.i][n];
                    H += J2*spins[j];
                }
            }

            this->mut.dS = -2*spins[mut.i] + 2.*spins[mut.i].dot(H)/std::pow(H.norm(),2) * H;
        }

        void generate_mutation() {
#ifdef CLUSTER_UPDATE
            SpinModel::generate_mutation(); 
#else
            mut.i = r() % V;
            mut_counter++;

            if (mut_counter == V) {
                mut_counter = 0;
                mut_mode++;
            }

            if (mut_mode < 10) {
                over_relaxation_mutation();
            } else if (mut_mode < 14) {
                metropolis_mutation();
            } else {
                metropolis_mutation();
                mut_mode = 0;
            }
#endif
        }

        const float onsite_func(const Eigen::Vector3f &S) {
            float E = 0;

            // Onsite interactions
            E -= B.dot(S);

            float phi = std::atan2(S[1], S[0]);
            float theta;
            if (S[2] > 1.0) { theta = PI; }
            else if (S[2] < -1.0) { theta = -PI; }
            else { theta = std::acos(S[2]); }

            E += K1*S[2]*S[2];
            //E += K2*pow(S[0]*S[0]+S[1]*S[1],2);
            E += K3*std::cos(6*phi)*std::pow(sin(theta), 6); // Sixfold magnetocrystalline field

            return E;
        }
};


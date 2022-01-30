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
        Eigen::Vector3d B; 

        Eigen::Matrix3d R;
        int mut_counter;
        int mut_mode;

        GaussianDist *dist_r;


        TrigonalModel(int N, int L, float J1, float J2, float K1, float K2, float K3,
                                    Eigen::Vector3d B) : SpinModel(1, N, N, L) {

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
            this->mut_mode = 0;
            this->mut_counter = 0; 

            std::function<float(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc = 
            [J1](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
                return -J1*S1.dot(S2);
            };

            Eigen::Vector3d v1; v1 << 1., 0., 0.;
            Eigen::Vector3d v2; v2 << 0.5, std::sqrt(3)/2., 0.;
            Eigen::Vector3d v3; v3 << 0.5, -std::sqrt(3)/2., 0.;
            this->add_bond(1,0,0,0, v1, bondfunc);
            this->add_bond(-1,0,0,0, -v1, bondfunc);
            this->add_bond(0,1,0,0, v2, bondfunc);
            this->add_bond(0,-1,0,0, -v2, bondfunc);
            this->add_bond(1,-1,0,0, v3, bondfunc);
            this->add_bond(-1,1,0,0, -v3, bondfunc);

            std::function<float(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc_inter = 
            [J2](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
                return J2*S1.dot(S2);
            };

            if ((std::abs(J2) > 1e-6) && (L > 1)) {
                Eigen::Vector3d v4; v4 << 0., 0., 1.;
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

        const Eigen::Vector3d molecular_field(int i) {
            float x = spins[i][0];
            float y = spins[i][1];
            float z = spins[i][2];
            Eigen::Vector3d H; H << K3*(6*pow(x, 5) - 60*pow(x, 3)*pow(y, 2) + 30*x*pow(y, 4)),
                                    K3*(-6*pow(y, 5) + 60*pow(x, 2)*pow(y, 3) - 30*pow(x, 4)*y),
                                    K1*pow(spins[i][2], 2);
            H += B;

            int j;
            for (int n = 0; n < 6; n++) {
                j = neighbors[i][n];
                H -= J1*spins[j];
            }
            for (int n = 6; n < 8; n++) {
                j = neighbors[i][n];
                H += J2*spins[j];
            }

            return H;
        }

#ifdef CLUSTER_UPDATE
        void generate_mutation() {
            cluster_update(); 
        }
#else
        void over_relaxation_mutation() {
            Eigen::Vector3d H = B;

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
        }
#endif

        const double onsite_func(const Eigen::Vector3d &S) {
            double E = 0;

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

        void rotate_spin(int i, Eigen::Vector3d v, float p) {
            float vx = v[0]/v.norm(); float vy = v[1]/v.norm(); float vz = v[2]/v.norm();
            Eigen::Matrix3d R; R << cos(p) + vx*vx*(1 - cos(p)), vx*vy*(1 - cos(p)) - vz*sin(p), vx*vz*(1 - cos(p)) + vy*sin(p),
                                    vy*vx*(1 - cos(p)) + vz*sin(p), cos(p) + vy*vy*(1 - cos(p)), vy*vz*(1 - cos(p)) - vx*sin(p),
                                    vz*vx*(1 - cos(p)) - vy*sin(p), vz*vy*(1 - cos(p)) + vx*sin(p), cos(p) + vz*vz*(1 - cos(p));
            spins[i] = R*spins[i];
        }

        void dynamic_step(float dt) {
            std::vector<Eigen::Vector3d> H = std::vector<Eigen::Vector3d>(V);
            int j;

            // Compute local molecular fields
            for (int i = 0; i < V; i++) {
                H[i] = molecular_field(i);
            }

            float dT;
            float Hm;

            // Precess around local molecular field
            for (int i = 0; i < V; i++) {
                Hm = H[i].norm();
                dT = Hm*dt;
                spins[i] = cos(dT)*spins[i] + sin(dT)*H[i].cross(spins[i])/Hm + (1 - cos(dT))*H[i].dot(spins[i])*H[i]/pow(Hm, 2);
            }
        }

};


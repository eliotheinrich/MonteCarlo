#include <iostream>
#include <stdlib.h>
#include <ctime>
#include "../MonteCarlo.cpp"
#include "../Spin2DModel.cpp"
#include "../Utility.cpp"


using namespace std;

Matrix3f generate_interaction_matrix(float J, Vector3f DM) {
    Matrix3f int_mat = Matrix3f::Constant(0);
    int_mat << -J, -DM[2], DM[1],
               DM[2], -J, -DM[0],
               -DM[1], DM[0], -J;
    return int_mat;
}

inline Matrix3f rotation_z(float t) {
    Matrix3f Rz = Matrix3f::Constant(0);
    Rz << cos(t), -sin(t), 0,
          sin(t), cos(t),  0,
          0,      0,       1;
    return Rz;
}

inline Vector3f z_comp(float t) {
    Vector3f v = Vector3f::Constant(0);
    v << 0, 0, cos(t);
    return v;
}

class HexagonalModel : public Spin2DModel<2> {
    public:
        float J1; float J2; float J3;
        float J4; float J5; float J6; float J7; float J8; float J9;
        float D1; float D2; float D3; float phi;
        float D4; float D5; float D6; float D7; float D8; float D9;
        float t4; float t5; float t6; float t7; float t8; float t9;
        float p4; float p5; float p6; float p7; float p8; float p9;
        float K1; float K2; float K3;
        float B;  float A;  float S;

        Matrix3f int_mat1;
        Matrix3f int_mat2;
        Matrix3f int_mat3;
        Matrix3f int_mat4;
        Matrix3f int_mat5;
        Matrix3f int_mat6;
        Matrix3f int_mat7;
        Matrix3f int_mat8;
        Matrix3f int_mat9;

        HexagonalModel(int N, float J1, float J2, float J3, 
                              float J4, float J5, float J6, float J7, float J8, float J9,
                              float D1, float D2, float D3, float phi, 
                              float D4, float t4, float p4, float D5, float t5, float p5, float D6, float t6, float p6,
                              float D7, float t7, float p7, float D8, float t8, float p8, float D9, float t9, float p9,
                              float K1, float K2, float K3, 
                              float B,  float A,  float S) : Spin2DModel<2>(N) {


            this->J1 = J1; this->J2 = J2; this->J3 = J3; 
            this->J4 = J4; this->J5 = J5; this->J6 = J6; this->J7 = J7; this->J8 = J8; this->J9 = J9;
            this->D1 = D1; this->D2 = D2; this->D3 = D3; this->phi = phi; 
            this->D4 = D4; this->D5 = D5; this->D6 = D6; this->D7 = D7; this->D8 = D8; this->D9 = D9;
            this->t4 = t4; this->t5 = t5; this->t6 = t6; this->t7 = t7; this->t8 = t8; this->t9 = t9;
            this->p4 = p4; this->p5 = p5; this->p6 = p6; this->p7 = p7; this->p8 = p8; this->p9 = p9;
            this->K1 = K1; this->K2 = K2; this->K3 = K3;
            this->B = B;   this->A = A;   this->S = S;

            Vector3f DMI1;
            DMI1 << D1*(-0.5)*cos(phi), D1*sqrt(3)/2*cos(phi), D1*sin(phi);
            Vector3f DMI2;
            DMI2 << D2*(-0.5)*cos(phi), -D2*sqrt(3)/2*cos(phi), D2*sin(phi);
            Vector3f DMI3;
            DMI3 << D3*cos(phi), 0, D3*sin(phi);


            Vector3f B1; B1 << sqrt(3), 0, 0;
            B1 = B1.normalized();
            Vector3f B2; B2 << sqrt(3)/2., 3/2., 0;
            B2 = B2.normalized();
            Vector3f B3; B3 << -sqrt(3)/2., 3/2., 0;
            B3 = B3.normalized();

            Vector3f DMI4 = D4*(rotation_z(p4)*B1*sin(t4) + z_comp(t4));
            Vector3f DMI5 = D5*(rotation_z(p5)*B2*sin(t5) - z_comp(t5));
            Vector3f DMI6 = D6*(rotation_z(p6)*B3*sin(t6) + z_comp(t6));

            Vector3f DMI7 = D7*(rotation_z(p7)*B1*sin(t7) - z_comp(t7));
            Vector3f DMI8 = D8*(rotation_z(p8)*B2*sin(t8) + z_comp(t8));
            Vector3f DMI9 = D9*(rotation_z(p9)*B3*sin(t9) - z_comp(t9));

            Vector3f g1; g1 << -1/sqrt(6), 1/sqrt(2), 1/sqrt(3);
            Vector3f g2; g2 << -1/sqrt(6), -1/sqrt(2), 1/sqrt(3);
            Vector3f g3; g3 << sqrt(2./3), 0., 1/sqrt(3);

            // NN interaction matrices
            this->int_mat1 = S*S*(generate_interaction_matrix(J1, DMI1) - K1 * g1 * g1.transpose());
            this->int_mat2 = S*S*(generate_interaction_matrix(J2, DMI2) - K2 * g2 * g2.transpose());
            this->int_mat3 = S*S*(generate_interaction_matrix(J3, DMI3) - K3 * g3 * g3.transpose());

            // NNN interaction matrices
            this->int_mat4 = S*S*generate_interaction_matrix(J4, DMI4);
            this->int_mat5 = S*S*generate_interaction_matrix(J5, DMI5);
            this->int_mat6 = S*S*generate_interaction_matrix(J6, DMI6);
            this->int_mat7 = S*S*generate_interaction_matrix(J7, DMI7);
            this->int_mat8 = S*S*generate_interaction_matrix(J8, DMI8);
            this->int_mat9 = S*S*generate_interaction_matrix(J9, DMI9);
        }

        const float onsite_energy(int i, int j, int s) {
            float E = 0.;

            E -= B*S*spins[i][j][s][2];
            E -= A*S*S*spins[i][j][s][2]*spins[i][j][s][2];

            return E;
        }

        const float bond_energy(int i, int j, int s) {
            float E = 0.;

            int low_i = mod(i-1, N);
            int low_j = mod(j-1, N);
            int high_i = (i+1)%N;
            int high_j = (j+1)%N;

            if (s == 0) {
                E += 0.5*float(spins[i][j][0].transpose() * (int_mat1 * spins[i][low_j][1]
                                                           + int_mat2 * spins[high_i][low_j][1]
                                                           + int_mat3 * spins[i][j][1]));

                E += 0.5*float(spins[i][j][0].transpose() * (int_mat4 * spins[high_i][j][0]
                                                           + int_mat4.transpose() * spins[low_i][j][0]
                                                           + int_mat5 * spins[i][high_j][0]
                                                           + int_mat5.transpose() * spins[i][low_j][0]
                                                           + int_mat6 * spins[low_i][high_j][0]
                                                           + int_mat6.transpose() * spins[high_i][low_j][0]));
            } else {
                E += 0.5*float(spins[i][j][1].transpose() * (int_mat1.transpose() * spins[i][high_j][0]
                                                           + int_mat2.transpose() * spins[low_i][high_j][0]
                                                           + int_mat3.transpose() * spins[i][j][0]));

                E += 0.5*float(spins[i][j][1].transpose() * (int_mat7 * spins[high_i][j][1]
                                                           + int_mat7.transpose() * spins[low_i][j][1]
                                                           + int_mat8 * spins[i][high_j][1]
                                                           + int_mat8.transpose() * spins[i][low_j][1]
                                                           + int_mat9 * spins[low_i][high_j][1]
                                                           + int_mat9.transpose() * spins[high_i][low_j][1]));

            }

            return E;
        }

};


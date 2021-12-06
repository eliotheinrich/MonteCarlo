#include "HexagonalModel_NEW.cpp"
#include "HexagonalModel.cpp"

#include<iostream>
#include<Eigen/Dense>


using namespace std;
using namespace Eigen;


int main() {
    int N = 20;
    float J1 = 0.2; float J2 = 0.4; float J3 = 0.3;
    float J4 = 0.1; float J5 = 0.6; float J6 = 0.1; 
    float J7 = 0.0; float J8 = 0.3; float J9 = 0.7;
    float D1 = 0.1; float D2 = 0.2; float D3 = 0.3; float phi = 0.5;
    float D4 = 1.; float t4 = .1; float p4 = .2; 
    float D5 = .6; float t5 = .3; float p5 = .4;
    float D6 = 0.5; float t6 = .8; float p6 = .9; 
    float D7 = 0.3; float t7 = 0.4; float p7 = -.5;
    float D8 = 0.2; float t8 = -.2; float p8 = -.2; 
    float D9 = 0.9; float t9 = -.5; float p9 = -.5;
    float K1 = 0.2; float K2 = 0.3; float K3 = 1.2;
    float B = 0.2;  float A = 0.22;  float S = 3.;

    srand(5);
    HexagonalModel model1 = HexagonalModel(N, J1, J2, J3, 
                                          J4, J5, J6, J7, J8, J9,
                                          D1, D2, D3, phi,
                                          D4, t4, p4, D5, t5, p5, D6, t6, p6,
                                          D7, t7, p7, D8, t8, p8, D9, t9, p9,
                                          K1, K2, K3, 
                                          B, A, S);
    srand(5);
    HexagonalModel_OLD model2 = HexagonalModel_OLD(N, J1, J2, J3, 
                                          J4, J5, J6, J7, J8, J9,
                                          D1, D2, D3, phi,
                                          D4, t4, p4, D5, t5, p5, D6, t6, p6,
                                          D7, t7, p7, D8, t8, p8, D9, t9, p9,
                                          K1, K2, K3, 
                                          B, A, S);

    model1.mut.dS << 0.2, 0.4, 0.3;
    model1.mut.i = 5;
    model1.mut.j = 2;
    model1.mut.s = 0;
    model2.mut.dS << 0.2, 0.4, 0.3;
    model2.mut.i = 5;
    model2.mut.j = 2;
    model2.mut.s = 0;

    float E1 = model1.energy();
    float dE = model1.energy_change();
    float E2 = model1.energy();


    float E1p = model2.energy();
    model2.spins[5][2][0] += model2.mut.dS;
    float E2p = model2.energy();


    cout << E1 << "\t\t" << E1p << endl;
    cout << E2 - E1 << endl;
    cout << E2p - E1p << endl;



}

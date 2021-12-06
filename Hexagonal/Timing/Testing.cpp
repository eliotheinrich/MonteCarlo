#include<iostream>
#include "../HexagonalModel.cpp"
#include "../../Utility.cpp"

void f1() {
    int N = 20;
    float J1 = 0.5; float J2 = 0.5; float J3 = 0.5;
    float J4 = 0.; float J5 = 0.; float J6 = 0.; 
    float J7 = 0.; float J8 = 0.; float J9 = 0.;
    float D1 = 0.; float D2 = 0.; float D3 = 0.; float phi = 0.;
    float D4 = 0.; float t4 = 0.; float p4 = 0.; 
    float D5 = 0.; float t5 = 0.; float p5 = 0.;
    float D6 = 0.; float t6 = 0.; float p6 = 0.; 
    float D7 = 0.; float t7 = 0.; float p7 = 0.;
    float D8 = 0.; float t8 = 0.; float p8 = 0.; 
    float D9 = 0.; float t9 = 0.; float p9 = 0.;
    float K1 = 0.; float K2 = 0.; float K3 = 0.;
    float B = 0.1;  float A = 0.2;  float S = 1.5;

    HexagonalModel model = HexagonalModel(N, J1, J2, J3, 
                                          J4, J5, J6, J7, J8, J9,
                                          D1, D2, D3, phi, 
                                          D4, t4, p4, D5, t5, p5, D6, t6, p6,
                                          D7, t7, p7, D8, t8, p8, D9, t9, p9,
                                          K1, K2, K3, 
                                          B, A, S);
    srand((unsigned)time( NULL ));
    MonteCarlo<HexagonalModel> m(&model);

    int nsteps = 100000;
    float T_max = .1;
    float T_min = 0.01;
    m.steps(nsteps, "trig", T_max, T_min); 
    cout << model.energy() << endl;
}

void f2() {
    int N = 20;
    float J1 = 0.5; float J2 = 0.5; float J3 = 0.5;
    float J4 = 0.; float J5 = 0.; float J6 = 0.; 
    float J7 = 0.; float J8 = 0.; float J9 = 0.;
    float D1 = 0.; float D2 = 0.; float D3 = 0.; float phi = 0.;
    float D4 = 0.; float t4 = 0.; float p4 = 0.; 
    float D5 = 0.; float t5 = 0.; float p5 = 0.;
    float D6 = 0.; float t6 = 0.; float p6 = 0.; 
    float D7 = 0.; float t7 = 0.; float p7 = 0.;
    float D8 = 0.; float t8 = 0.; float p8 = 0.; 
    float D9 = 0.; float t9 = 0.; float p9 = 0.;
    float K1 = 0.; float K2 = 0.; float K3 = 0.;
    float B = 0.1;  float A = 0.2;  float S = 1.5;

    HexagonalModel model = HexagonalModel(N, J1, J2, J3, 
                                          J4, J5, J6, J7, J8, J9,
                                          D1, D2, D3, phi, 
                                          D4, t4, p4, D5, t5, p5, D6, t6, p6,
                                          D7, t7, p7, D8, t8, p8, D9, t9, p9,
                                          K1, K2, K3, 
                                          B, A, S);
    srand((unsigned)time( NULL ));
    MonteCarlo<HexagonalModel> m(&model);

    int nsteps = 1000;
    float T_max = 0.1;
    float T_min = 0.01;
    for (int i = 0; i < 100; i++) {
        m.steps(nsteps, "", T_max, T_min); 
    }
    cout << model.energy() << endl;
}


int main() {
    cout << "Method 1: " << time_code(f1) << endl;
    cout << "Method 2: " << time_code(f2) << endl;

}






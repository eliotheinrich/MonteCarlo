#include<Eigen/Dense>
#include<iostream>
#include "TrigonalModel.cpp"

using namespace std;
using namespace Eigen;

int main() {
    int N = 8;
    int L = 1;

    float J1 = 1.;
    float J2 = 0.;
    float K1 = 0.2;
    float K2 = 0.;
    float K3 = 0.;

    Vector3f a1; a1 << 1, 0, 0;
    Vector3f a2; a2 << 0.5, 0.866025, 0;
    Vector3f ab = (a1 + a2).normalized();


    float Bm = 0.1;
    Vector3f B; B << 1, 0, 0;
    B = B*Bm;
    float S = 1.;

    const int layers = 2;

    srand((unsigned)time( NULL ));
    TrigonalModel model = TrigonalModel(N, L, J1, J2, K1, K2, K3, B, S);


    int MCSteps = N*N*L;
    int nsteps = 1000*MCSteps;
    float T = 0.1;
    run_MC<MagnetizationLogItem>(&model, nsteps, "const", T, T, true);

    model.save_spins("SpinTexture.txt");
    cout << model.energy() << endl;


}



#include<Eigen/Dense>
#include<iostream>
#include "SquareModel.cpp"

using namespace std;
using namespace Eigen;

int main() {
    int N = 8;
    int L = 1;

    float J = 1.;
    float A = 0.;

    srand((unsigned)time( NULL ));
    SquareModel model = SquareModel(N, L, J, A);


    int MCSteps = N*N*L;
    int nsteps = 1000*MCSteps;
    float T = 0.1;
    run_MC<MagnetizationLogItem>(&model, nsteps, "const", T, T, true);

    model.save_spins("SpinTexture.txt");
    cout << model.energy() << endl;


}



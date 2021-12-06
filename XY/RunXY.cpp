#include "SquareXYModel.cpp"
#include <iostream>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 8;
    const int L = 1;
    const float J = 1.;
    const float B = 0.;
    const float Bp = 0.;
    const float T = 0.01;


    const int MCStep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);

    run_MC<MagnetizationLogItem>(model, 500*MCStep, "trig", T, T, true);
    cout << model->energy() << endl;

    model->save_spins("Spins.txt");
    

}

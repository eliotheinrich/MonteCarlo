#include "SquareIsingModel.cpp"
#include <iostream>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 20;
    const int L = 1;
    const float J = 1.;
    const float B = 0.;
    const float T = 0.1;


    const int MCStep = N*N*L;

    SquareIsingModel *model = new SquareIsingModel(N, L, J, B);

    cout << model->get_magnetization() << endl;
    run_MC<MagnetizationLogItem>(model, 50000*MCStep, "trig", T, T, true);
    cout << model->get_magnetization() << endl;

    model->save_spins("Spins.txt");
    

}

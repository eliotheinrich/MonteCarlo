#include "EasyPlaneHeis.cpp"
#include <iostream>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 8;
    const int L = 1;
    const float J = 1.;
    const float K = 0.3;
    const float T = 0.1;


    const int MCStep = N*N*L;

    EasyPlaneHeis *model = new EasyPlaneHeis(N, L, J, K);

    MonteCarlo<EasyPlaneHeis> *m = new MonteCarlo<EasyPlaneHeis>(model);

    function<vector<float>(EasyPlaneHeis*)> logfunc = [](EasyPlaneHeis* model) { return MagnetizationLogItem(model); };
    vector<vector<float>> log = run_MC(m, 500*MCStep, T, logfunc, 100);
    write_to_file("Log.txt", log);


    cout << "Energy = " << model->energy() << endl;

    model->save_spins("Spins.txt");
    

}

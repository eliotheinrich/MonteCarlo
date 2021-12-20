#include "TrigonalModel.cpp"
#include <iostream>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 8;
    const int L = 1;
    const float J1 = 1.;
    const float J2 = 0.01;
    const float K1 = 0.5;
    const float K2 = 0.;
    const float K3 = 0.03;
    Vector3f B; B << 0., 0., 0.;
    const float T = 0.1;


    const int MCStep = N*N*L;

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
    MonteCarlo<TrigonalModel> *m = new MonteCarlo<TrigonalModel>(model);

    function<vector<float>(TrigonalModel*)> logfunc = [](TrigonalModel* model) { return MagnetizationLogItem(model); };
    vector<vector<float>> log = run_MC(m, 500*MCStep, T, logfunc, 100);
    write_to_file("Log.txt", log);

    cout << "Energy = " << model->energy() << endl;

}

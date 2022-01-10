#include "TrigonalModel.cpp"
#include <iostream>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 8;
    const int L = 1;
    const float J1 = 1.;
    const float J2 = 0.0;
    const float K1 = 0.;
    const float K2 = 0.;
    const float K3 = 0.0;
    Vector3f B; B << 0., 0., 0.;
    const float T = .4;


    const int MCStep = N*N*L;

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
    MonteCarlo<TrigonalModel> *m = new MonteCarlo<TrigonalModel>(model);

    int num_samples = 1000;
    vector<vector<float>> log = vector<vector<float>>(num_samples, vector<float>(2));
    for (int i = 0; i < num_samples; i++) {
        m->steps(MCStep, T);
        log[i][0] = model->get_magnetization().norm();
        log[i][1] = model->energy();
    }

    ofstream output("Log.txt");
    for (int i = 0; i < log.size(); i++) {
        output << log[i][0] << "\t" << log[i][1] << "\n";
    }

    model->save_spins("Spins.txt");

    cout << "Energy = " << model->energy() << endl;

}

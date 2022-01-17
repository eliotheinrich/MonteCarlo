#include "SquareClockModel.cpp"
#include <iostream>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 8;
    const int L = 1;
    const float J = 1.;
    const float B = 0.;
    const float Bp = 0.;
    const float T = 0.1;


    const int MCStep = 1;//N*N*L;
    const int q = 6;

    SquareClockModel<q> *model = new SquareClockModel<q>(N, L, J);
    MonteCarlo<SquareClockModel<q>> *m = new MonteCarlo<SquareClockModel<q>>(model);

    int num_samples = 1000;
    vector<vector<float>> log = vector<vector<float>>(num_samples, vector<float>(2));
    for (int i = 0; i < num_samples; i++) {
        m->steps(MCStep, T);
        log[i][0] = model->get_magnetization();
        log[i][1] = model->energy();
    }

    ofstream output("Log.txt");
    for (int i = 0; i < log.size(); i++) {
        output << log[i][0] << "\t" << log[i][1] << "\n";
    }

    model->save_spins("Spins.txt");
    //model->cluster_update();
    //cout << model->s.size() << endl;
    //model->save_spins("Spins2.txt");

    cout << "Energy = " << model->energy() << endl;

}

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

    MonteCarlo<SquareXYModel> *m = new MonteCarlo<SquareXYModel>(model);

    function<vector<float>(SquareXYModel*)> logfunc = [](SquareXYModel* model) { return MagnetizationLogItem(model); };
    vector<vector<float>> log = run_MC(m, 500*MCStep, T, logfunc, 100);
    write_to_file("Log.txt", log);


    function<float(SquareXYModel*)> stiffness = [T](SquareXYModel *m) { return m->spin_stiffness(T); };
    float E = m->expectation(stiffness, T, 100, 10);
    cout << "Energy = " << E << endl;

    cout << model->energy() << endl;

    model->save_spins("Spins.txt");
    

}

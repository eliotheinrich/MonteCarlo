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
    SquareModel *model = new SquareModel(N, L, J, A);
    MonteCarlo<SquareModel> *m = new MonteCarlo<SquareModel>(model);


    int MCStep = N*N*L;
    int nsteps = 1000*MCStep;
    float T = 0.1;

    function<vector<float>(SquareModel*)> logfunc = [](SquareModel* model) { return MagnetizationLogItem(model); };
    vector<vector<float>> log = run_MC(m, 500*MCStep, T, logfunc, 100);
    write_to_file("Log.txt", log);

    model->save_spins("SpinTexture.txt");
    cout << model->energy() << endl;


}



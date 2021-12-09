#include "TrigonalModel.cpp"
#include <iostream>
#include <chrono>

using namespace std;
using namespace Eigen;


int main() {
    const int N = 20;
    const int L = 4;

    const float J1 = 1.5;
    const float J2 = 0.4;
    const float K1 = 0.8;
    const float K2 = 0.1;
    const float K3 = 0.1;

    const float T = 0.155;

    Vector3f Bhat1; Bhat1 << 0.866025, 0.5, 0;
    Vector3f Bhat2; Bhat2 << 1., 0., 0.;
    float Bm;
  
    const int MCstep = N*N*L;

    srand((unsigned)time( NULL ));
    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, Bm*Bhat1);
    MonteCarlo<TrigonalModel> *m = new MonteCarlo<TrigonalModel>(model);


    int nsteps = 10000*MCstep;


    auto start = chrono::high_resolution_clock::now();
    run_MC(m, nsteps, "trig", 3*J1, T);
    auto stop = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int microseconds = duration.count();

    cout << to_string(nsteps) << " steps took " << to_string(float(microseconds)/1000000.) << "s." << endl;
    cout << to_string(float(nsteps)/(float(microseconds)/1000000)) << " steps/second." << endl;



}


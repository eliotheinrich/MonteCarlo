#include "SquareXYModel.cpp"
#include<iostream>
#include<chrono>

using namespace std;
using namespace Eigen;


int main() {

    srand((unsigned)time( NULL ));

    const int N = 20;
    const int L = 4;

    const float J = 1.;
    const float B = 0.;
    const float Bp = 0.;

    const float T = 0.155;
  
    const int MCstep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);
    model->random_selection = false;


    int nsteps = 1000*MCstep;


    auto start = chrono::high_resolution_clock::now();
    run_MC(model, nsteps, "trig", T, T);
    auto stop = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int microseconds = duration.count();

    cout << to_string(nsteps) << " steps took " << to_string(float(microseconds)/1000000.) << "s." << endl;
    cout << to_string(float(nsteps)/(float(microseconds)/1000000)) << " steps/second." << endl;

//    cout << "over-relax steps: " << model->f1 << endl;
//    cout << "standard steps: " << model->f2 << endl;



}


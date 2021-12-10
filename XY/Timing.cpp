#include "SquareXYModel.cpp"
#include <iostream>
#include <chrono>
#include <ctpl_stl.h>

using namespace std;
using namespace Eigen;


int main(int charc, char* argv[]) {

    srand((unsigned)time( NULL ));

    const int N = 20;
    const int L = 4;

    const float J = 1.;
    const float B = 0.;
    const float Bp = 0.;

    const float T = 0.155;
  
    const int MCstep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);
    MonteCarlo<SquareXYModel> *m = new MonteCarlo<SquareXYModel>(model);


    //int num_trials = stoi(argv[1]);
    int num_trials = 2;
    cout << "Num_trials = " << num_trials << endl;
    int nsteps = 1000*MCstep;

    //ctpl::thread_pool threads(2);



    auto func = [&m, nsteps, J, T](int id) { 
        run_MC(m, nsteps, "trig", J, T); 
    };

    auto start = chrono::high_resolution_clock::now();

    //vector<future<void>> results(num_trials);

    //for (int i = 0; i < num_trials; i++) {
    //    results[i] = threads.push(func);
    //}

    //for (int i = 0; i < num_trials; i++) {
    //    results[i].get();
    //}
    thread p1(func, 1);
    thread p2(func, 1);


    p1.join();
    p2.join();

    auto stop = chrono::high_resolution_clock::now();
    

    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int microseconds = duration.count();

    cout << to_string(nsteps*num_trials) << " steps took " << to_string(float(microseconds)/1000000.) << "s." << endl;
    cout << to_string(float(nsteps)/(float(microseconds)/1000000)) << " steps/second." << endl;

}


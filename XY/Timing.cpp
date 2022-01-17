#include "TrigonalXYModel.cpp"
#include <iostream>
#include <chrono>
#include <ctpl.h>

int main(int charc, char* argv[]) {

    srand((unsigned)time( NULL ));

    const int N = 20;
    const int L = 4;

    const float J = 1.;
    const float A = 0.2;

    const float T = 0.155;
  
    const int MCstep = N*N*L;

    TrigonalXYModel *model = new TrigonalXYModel(N, L, J, A);
    MonteCarlo<TrigonalXYModel> *m = new MonteCarlo<TrigonalXYModel>(model);


    //int num_trials = stoi(argv[1]);
    int nsteps = 1000*MCstep;

    auto start = std::chrono::high_resolution_clock::now();

    m->steps(nsteps, T);

    auto stop = std::chrono::high_resolution_clock::now();
    

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int microseconds = duration.count();

    std::cout << to_string(nsteps) << " steps took " << to_string(float(microseconds)/1000000.) << "s." << std::endl;
    std::cout << to_string(float(nsteps)/(float(microseconds)/1000000)) << " steps/second." << std::endl;

}


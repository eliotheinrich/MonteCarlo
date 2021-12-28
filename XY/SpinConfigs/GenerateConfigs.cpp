#include "../TrigonalXYModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    string foldername = argv[1];
    int N = stoi(argv[2]);
    const int L = 1;
    const float J = 1.;
    const float A = 0.3;

    const int MCStep = N*N*L;
    cout << "N = " << N << endl;

    TrigonalXYModel *model = new TrigonalXYModel(N, L, J, A);

    const float Tmax = 3.;
    const float Tmin = 0.1;
    unsigned long long int resolution = stoi(argv[3]);

    vector<float> T(resolution);

    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    unsigned long long int equilibration_steps = stoi(argv[4])*MCStep;
    unsigned long long int num_samples = stoi(argv[5]);
    unsigned long long int steps_per_sample = stoi(argv[6])*MCStep;

    int num_threads = stoi(argv[7]);

    auto start = chrono::high_resolution_clock::now();

    generate_spin_configs(model, T, equilibration_steps, num_samples, steps_per_sample, num_threads, foldername);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    cout << "Completion time: " << seconds/60. << " minutes." << endl;


}


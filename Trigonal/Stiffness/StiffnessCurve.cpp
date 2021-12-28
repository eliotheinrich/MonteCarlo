#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

#define BOLTZMANN_CONSTANT 0.08617

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    string filename = argv[1];
    int N = stoi(argv[2]);
    int num_threads = stoi(argv[3]);

    cout << "num_threads = " << num_threads << endl;
    cout << "N = " << N << endl;

    const int L = N/2;
    const float J1 = 0.64;
    const float J2 = 0.;
    const float K1 = 2.56;
    const float K2 = 0.;
    const float K3 = 0.15;
    Vector3f B; B << 0., 0., 0.3;

    const int MCStep = N*N*L;

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);

    const float Tmax = 3.;
    const float Tmin = 0.1;
    unsigned long long int resolution = 30;

    vector<float> T(resolution);
    ofstream output(filename);

    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    unsigned long long int num_runs = 5;
    unsigned long long int steps_per_run = 5000*MCStep;

    unsigned long long int num_samples = 1000;
    unsigned long long int steps_per_sample = 30*MCStep;

    auto start = chrono::high_resolution_clock::now();

    stiffness_run(model, &T, num_runs, steps_per_run, num_samples, steps_per_sample, num_threads, filename);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    cout << "Completion time: " << seconds/60. << " minutes." << endl;


}


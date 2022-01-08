#include "../TrigonalXYModel.cpp"
#include "../SquareXYModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

using namespace std;
using namespace Eigen;

vector<double> twist_sampler(SquareXYModel *model) {
    return model->twist_stiffness();
}

int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    string filename = argv[1];
    int N = stoi(argv[2]);
    cout << "N = " << N << endl;
    int num_threads = 4;

    const int L = 1;
    const float J = 1.;
    const float A = 0.3;

    const int MCStep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, 0., 0.);

    unsigned long long int resolution = 30;
    unsigned long long int num_runs = 5;
    unsigned long long int steps_per_run = 5000*MCStep;

    unsigned long long int num_samples = 100;
    unsigned long long int steps_per_sample = 100*MCStep;

    vector<float> T(resolution);

    const float Tmax = 3.;
    const float Tmin = 0.1;
    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    auto start = chrono::high_resolution_clock::now();

    auto data = sample_r(twist_sampler, model, &T, 1, steps_per_run, num_samples, steps_per_sample, num_threads);
    auto stats = summary_statistics(&data);
    write_data(&stats, &T, filename);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    cout << "Completion time: " << seconds/60. << " minutes." << endl;


}


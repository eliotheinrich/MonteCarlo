#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

#define BOLTZMANN_CONSTANT 0.08617

using namespace std;
using namespace Eigen;

vector<double> twisting_sampler(TrigonalModel *model) {
    return model->twist_stiffness();
}

int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    string filename = argv[1];
    int num_threads = stoi(argv[2]);
    int N = stoi(argv[3]);
    int L = 1;

    cout << "N = " << N << endl;

    // Base case
    float J1 = stof(argv[4]);
    float J2 = stof(argv[5]);
    float K1 = stof(argv[6]);
    float K2 = stof(argv[7]);
    float K3 = stof(argv[8]);
    Vector3f B; B << stof(argv[9]), stof(argv[10]), stof(argv[11]);

    const int MCStep = N*N*L;

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);

    const float Tmax = 3.;
    const float Tmin = 0.1;
    unsigned int resolution = stoi(argv[12]);

    vector<float> T(resolution);

    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    unsigned long long int steps_per_run = stoi(argv[13])*MCStep;

    unsigned int num_samples = stoi(argv[14]);
    unsigned long long int steps_per_sample = stoi(argv[15])*MCStep;

    int num_runs = stoi(argv[16]);

    unsigned long long nsteps = resolution*num_runs*(steps_per_run + num_samples*steps_per_sample);
    cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << endl;

    auto start = chrono::high_resolution_clock::now();

    auto data = sample_pt(twisting_sampler, model, &T, steps_per_run, num_samples, steps_per_sample, num_threads);
    auto stats = summary_statistics(&data);
    write_data(&data, &T, filename);

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    cout << "Completion time: " << seconds/60. << " minutes." << endl;


}


#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

using namespace std;
using namespace Eigen;


vector<float> Cij_sampler(TrigonalModel *model) {
    vector<float> Cij = vector<float>(model->V, 0.);
    vector<float> Cij_avg = vector<float>(model->V, 0.);

    int i;
    int j;
    for (int i = 0; i < model->V; i++) {
        Cij = model->skyrmion_correlation_function(i);
        for (int j = 0; j < model->V; j++) {
            Cij_avg[j] += Cij[j];
        }
    }

    for (int i = 0; i < model->V; i++) {
        Cij_avg[i] = Cij_avg[i]/model->V;
    }

    return Cij_avg;
}

int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    string filename = argv[1];
    int num_threads = stoi(argv[2]);
    int N = stoi(argv[3]);
    int L = stoi(argv[4]);

    cout << "N = " << N << endl;
    cout << "L = " << N << endl;

    // Base case
    float J1 = stof(argv[5]);
    float J2 = stof(argv[6]);
    float K1 = stof(argv[7]);
    float K2 = stof(argv[8]);
    float K3 = stof(argv[9]);
    Vector3f B; B << stof(argv[10]), stof(argv[11]), stof(argv[12]);

    const int MCStep = N*N*L;

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);

    const float Tmax = 3.;
    const float Tmin = 0.1;
    unsigned int resolution = stoi(argv[13]);

    vector<float> T(resolution);

    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    unsigned long long int steps_per_run = stoi(argv[14])*MCStep;

    unsigned int num_samples = stoi(argv[15]);
    unsigned long long int steps_per_sample = stoi(argv[16])*MCStep;

    int num_runs = stoi(argv[17]);

    unsigned long long nsteps = resolution*num_runs*(steps_per_run + num_samples*steps_per_sample);
    cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << endl;

    auto start = chrono::high_resolution_clock::now();

    auto data = sample_r(Cij_sampler, model, &T, num_runs, steps_per_run, num_samples, steps_per_sample, num_threads);
    auto stats = summary_statistics(&data);
    string header = to_string(model->N1) + "\t" + to_string(model->N2) + "\t" + to_string(model->N3);
    write_data(&stats, &T, filename, header);


    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    cout << "Completion time: " << seconds/60. << " minutes." << endl;


}


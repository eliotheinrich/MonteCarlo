#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

using namespace std;
using namespace Eigen;

vector<float> magnetization_sampler(TrigonalModel *model) {
    return vector<float>{
                static_cast<float>(model->get_magnetization().norm()),
                static_cast<float>(model->energy()/model->V)
           };
}


int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    // Folder to store data
    string foldername = argv[1];

    // Load params config file
    ifstream paramfile(argv[2]);

    //float S = 3.5;
    float S = 1.;

    // Load paramfile line by line
    string line;
    getline(paramfile, line);
    getline(paramfile, line);
    vector<string> paramss = split(&line, ",");

    int N = stoi(paramss[0]);
    int L = stoi(paramss[1]);
    const int MCStep = N*N*L;

    float J1 = S*S*stof(paramss[2]);
    float J2 = S*S*stof(paramss[3]);
    float K1 = S*stof(paramss[4]);
    float K2 = S*stof(paramss[5]);
    float K3 = S*stof(paramss[6]);
    float Bx = stof(paramss[7]);
    float By = stof(paramss[8]);
    float Bz = stof(paramss[9]);
    Vector3f B; B << Bx, By, Bz;

    getline(paramfile, line);
    getline(paramfile, line);
    getline(paramfile, line);
    paramss = split(&line, ",");

    int resolution = stoi(paramss[0]);
    int steps_per_run = stoi(paramss[1]);
    int num_samples = stoi(paramss[2]);
    int steps_per_sample = stoi(paramss[3]);
    int num_runs = stoi(paramss[4]);

    //float T_max = 60*BOLTZMANN_CONSTANT; // In Kelvin
    //float T_min = 0.1*BOLTZMANN_CONSTANT;
    float T_max = 2.;
    float T_min = 0.05;
    vector<float> T(resolution);
    for (int i = 0; i < resolution; i++) {
        T[i] = T_max*i/resolution + T_min*(resolution - i)/resolution;
    }

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);

    int num_threads = stoi(argv[4]);
    unsigned long long int nsteps = (long long) resolution*(steps_per_run + num_samples*steps_per_sample);

    cout << "Number steps: " << nsteps << endl;
    cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << endl;

    auto start = chrono::high_resolution_clock::now();

    auto data = sample_r(magnetization_sampler, model, &T, 5, steps_per_run, num_samples, steps_per_sample, num_threads);
    auto stats = summary_statistics(&data);
    write_data(&stats, &T, foldername + "/MagnetizationCurve.txt");

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    cout << "Completion time: " << seconds/60. << " minutes." << endl;

}


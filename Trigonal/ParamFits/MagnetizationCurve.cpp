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

    // Magnetic field strength
    float Bm = S*stof(argv[3]);

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

    getline(paramfile, line);
    getline(paramfile, line);
    getline(paramfile, line);
    paramss = split(&line, ",");

    int resolution = stoi(paramss[0]);
    int steps_per_run = stoi(paramss[1])*MCStep;
    int num_samples = stoi(paramss[2]);
    int steps_per_sample = stoi(paramss[3])*MCStep;
    int num_runs = stoi(paramss[4]);

    //float T_max = 60*BOLTZMANN_CONSTANT; // In Kelvin
    //float T_min = 0.1*BOLTZMANN_CONSTANT;
    float T_max = 4.;
    float T_min = 0.05;
    vector<float> T(resolution);
    for (int i = 0; i < resolution; i++) {
        T[i] = T_max*i/resolution + T_min*(resolution - i)/resolution;
    }

    Vector3f Bhat1; Bhat1 << 0., 1., 0.;
    Vector3f Bhat2; Bhat2 << 1., 0., 0.;
    Vector3f Bhat3; Bhat3 << 0., 0., 1.;
    Vector3f B1 = Bm*Bhat1;
    Vector3f B2 = Bm*Bhat2; 
    Vector3f B3 = Bm*Bhat3;

    TrigonalModel *model1 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B1);
    TrigonalModel *model2 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B2);
    TrigonalModel *model3 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B3);

    int num_threads = stoi(argv[4]);
    unsigned long long int nsteps = (long long) 3*resolution*(steps_per_run + num_samples*steps_per_sample);

    cout << "Number steps: " << nsteps << endl;
    cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << endl;

    auto start = chrono::high_resolution_clock::now();

    auto data1 = sample_pt(magnetization_sampler, model1, &T, num_runs, steps_per_run, num_samples, steps_per_sample, num_threads);
    auto stats1 = summary_statistics(&data1)
    write_data(&stats1, &T, foldername + "/MagnetizationCurve1.txt");

    auto data2 = sample_pt(magnetization_sampler, model2, &T, num_runs, steps_per_run, num_samples, steps_per_sample, num_threads)
    auto stats2 = summary_statistics(&data2)
    write_data(&stats2, &T, foldername + "/MagnetizationCurve2.txt");

    auto data3 = sample_pt(magnetization_sampler, model3, &T, num_runs, steps_per_run, num_samples, steps_per_sample, num_threads)
    auto stats3 = summary_statistics(&data3)
    write_data(&stats3, &T, foldername + "/MagnetizationCurve3.txt");

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    cout << "Completion time: " << seconds/60. << " minutes." << endl;


}

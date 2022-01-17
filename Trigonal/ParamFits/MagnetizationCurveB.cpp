#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

std::vector<float> magnetization_sampler(TrigonalModel *model) {
    return std::vector<float>{
                static_cast<float>(model->get_magnetization().norm()),
                static_cast<float>(model->energy()/model->V)
           };
}


int main(int argc, char* argv[]) {
    bool cluster = true;
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
    std::vector<string> paramss = split(&line, ",");

    int N = std::stoi(paramss[0]);
    int L = std::stoi(paramss[1]);
    const int MCStep = N*N*L;

    float J1 = S*S*std::stof(paramss[2]);
    float J2 = S*S*std::stof(paramss[3]);
    float K1 = S*std::stof(paramss[4]);
    float K2 = S*std::stof(paramss[5]);
    float K3 = S*std::stof(paramss[6]);
    float Bx = std::stof(paramss[7]);
    float By = std::stof(paramss[8]);
    float Bz = std::stof(paramss[9]);
    Eigen::Vector3f B; B << Bx, By, Bz;

    getline(paramfile, line);
    getline(paramfile, line);
    getline(paramfile, line);
    paramss = split(&line, ",");

    int resolution = std::stoi(paramss[0]);
    int steps_per_run = std::stoi(paramss[1]);
    if (!cluster) { steps_per_run *= MCStep; }
    int num_samples = std::stoi(paramss[2]);
    int steps_per_sample = std::stoi(paramss[3]);
    if (!cluster) { steps_per_sample *= MCStep; }
    int num_runs = std::stoi(paramss[4]);

    //float T_max = 60*BOLTZMANN_CONSTANT; // In Kelvin
    //float T_min = 0.1*BOLTZMANN_CONSTANT;
    float T_max = 2.;
    float T_min = 0.05;
    std::vector<float> T(resolution);
    for (int i = 0; i < resolution; i++) {
        T[i] = T_max*i/resolution + T_min*(resolution - i)/resolution;
    }

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
    model->cluster = cluster;


    int num_threads = std::stoi(argv[4]);
    unsigned long long int nsteps = (long long) resolution*(steps_per_run + num_samples*steps_per_sample);

    std::cout << "Number steps: " << nsteps << std::endl;
    std::cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto data = sample_r(magnetization_sampler, model, &T, 5, steps_per_run, num_samples, steps_per_sample, num_threads);
    auto stats = summary_statistics(&data);
    write_data(&stats, &T, foldername + "/MagnetizationCurve.txt");

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    std::cout << "Completion time: " << seconds/60. << " minutes." << std::endl;

}


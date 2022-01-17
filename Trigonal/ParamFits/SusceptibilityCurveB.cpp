#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

std::vector<float> susceptibility_sampler(TrigonalModel *model) {
    return std::vector<float>{
                static_cast<float>(model->get_magnetization().norm()/pow(model->B.norm(), 2)),
                static_cast<float>(model->energy()/model->V)
           };
}


int main(int argc, char* argv[]) {
    std::srand((unsigned)std::time( NULL ));

    // Folder to store data
    std::string foldername = argv[1];

    // Load params config file
    std::ifstream paramfile(argv[2]);

    //float S = 3.5;
    float S = 1.;

    // Load paramfile line by line
    std::string line;
    getline(paramfile, line);
    getline(paramfile, line);
    std::vector<std::string> paramss = split(&line, ",");

    int N = std::stoi(paramss[0]);
    int L = std::stoi(paramss[1]);
    const int MCStep = N*N*L;

    float J1 = S*S*stof(paramss[2]);
    float J2 = S*S*stof(paramss[3]);
    float K1 = S*stof(paramss[4]);
    float K2 = S*stof(paramss[5]);
    float K3 = S*stof(paramss[6]);
    float Bx = stof(paramss[7]);
    float By = stof(paramss[8]);
    float Bz = stof(paramss[9]);
    Eigen::Vector3f B; B << Bx, By, Bz;

    getline(paramfile, line);
    getline(paramfile, line);
    getline(paramfile, line);
    paramss = split(&line, ",");

    int resolution = std::stoi(paramss[0]);
    int steps_per_run = std::stoi(paramss[1])*MCStep;
    int num_samples = std::stoi(paramss[2]);
    int steps_per_sample = std::stoi(paramss[3])*MCStep;
    int num_runs = std::stoi(paramss[4]);

    //float T_max = 60*BOLTZMANN_CONSTANT; // In Kelvin
    //float T_min = 0.1*BOLTZMANN_CONSTANT;
    float T_max = 4.;
    float T_min = 0.05;
    std::vector<float> T(resolution);
    for (int i = 0; i < resolution; i++) {
        T[i] = T_max*i/resolution + T_min*(resolution - i)/resolution;
    }

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);

    int num_threads = std::stoi(argv[4]);
    unsigned long long int nsteps = (long long) resolution*(steps_per_run + num_samples*steps_per_sample);

    std::cout << "Number steps: " << nsteps << std::endl;
    std::cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << std::endl;

    auto start = std::chrono::high_resolution_clock::now();


    auto data = sample_pt(susceptibility_sampler, model, &T, steps_per_run, num_samples, steps_per_sample, num_threads,
    write_data(&data, &T, foldername + "/SusceptibilityCurve.txt");

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    std::cout << "Completion time: " << seconds/60. << " minutes." << std::endl;

}


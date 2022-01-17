#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

#define BOLTZMANN_CONSTANT 0.08617

std::vector<double> twisting_sampler(TrigonalModel *model) {
    return model->twist_stiffness();
}

int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    std::string filename = argv[1];
    int num_threads = std::stoi(argv[2]);
    int N = std::stoi(argv[3]);
    int L = std::stoi(argv[4]);

    std::cout << "N = " << N << std::endl;
    std::cout << "L = " << L << std::endl;

    // Base case
    float J1 = std::stof(argv[5]);
    float J2 = std::stof(argv[6]);
    float K1 = std::stof(argv[7]);
    float K2 = std::stof(argv[8]);
    float K3 = std::stof(argv[9]);
    Eigen::Vector3f B; B << std::stof(argv[10]), std::stof(argv[11]), std::stof(argv[12]);

    const int MCStep = N*N*L;

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
    model->cluster = false;

    const float Tmax = 2.;
    const float Tmin = 0.1;
    unsigned int resolution = std::stoi(argv[13]);

    std::vector<float> T(resolution);

    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    unsigned long long steps_per_run = std::stoi(argv[14])*MCStep;

    unsigned int num_samples = std::stoi(argv[15]);
    unsigned long long steps_per_sample = std::stoi(argv[16])*MCStep;

    int num_runs = std::stoi(argv[17]);

    unsigned long long nsteps = resolution*num_runs*(steps_per_run + num_samples*steps_per_sample);
    std::cout << "nsteps: " << nsteps << std::endl;
    std::cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto data = sample_pt(twisting_sampler, model, &T, steps_per_run, num_samples, steps_per_sample, num_threads);
//    auto stats = summary_statistics(&data);
    std::string header = std::to_string(num_samples) + "\t" + std::to_string(N) + "\t" + std::to_string(L);
    write_samples(&data, &T, filename, header);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    std::cout << "Completion time: " << seconds/60. << " minutes." << std::endl;


}


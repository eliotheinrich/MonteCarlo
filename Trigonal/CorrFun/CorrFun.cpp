#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

std::vector<float> Cij_sampler(TrigonalModel *model) {
    std::vector<float> Cij = std::vector<float>(model->V, 0.);
    std::vector<float> Cij_avg = std::vector<float>(model->V, 0.);

    int i;
    int j;
    for (int i = 0; i < model->V; i++) {
        Cij = model->full_correlation_function(i);
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
    std::srand((unsigned)std::time( NULL ));

    std::string filename = argv[1];
    int num_threads = std::stoi(argv[2]);
    int N = std::stoi(argv[3]);
    int L = 1;

    std::cout << "N = " << N << std::endl;

    // Base case
    float J1 = std::stof(argv[4]);
    float J2 = std::stof(argv[5]);
    float K1 = std::stof(argv[6]);
    float K2 = std::stof(argv[7]);
    float K3 = std::stof(argv[8]);
    Eigen::Vector3f B; B << std::stof(argv[9]), stof(argv[10]), stof(argv[11]);

    const int MCStep = N*N*L;

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);

    const float Tmax = 3.;
    const float Tmin = 0.1;

    unsigned int resolution = std::stoi(argv[12]);

    std::vector<float> T(resolution);

    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    unsigned long long int steps_per_run = std::stoi(argv[13])*MCStep;

    unsigned int num_samples = std::stoi(argv[14]);
    unsigned long long int steps_per_sample = std::stoi(argv[15])*MCStep;

    int num_runs = std::stoi(argv[16]);

    unsigned long long nsteps = resolution*num_runs*(steps_per_run + num_samples*steps_per_sample);
    std::cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto data = sample_r(Cij_sampler, model, &T, num_runs, steps_per_run, num_samples, steps_per_sample, num_threads);
    auto stats = summary_statistics(&data);
    std::string header = std::to_string(model->N1) + "\t" + std::to_string(model->N2) + "\t" + std::to_string(model->N3);
    write_data(&stats, &T, filename, header);


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    std::cout << "Completion time: " << seconds/60. << " minutes." << std::endl;


}


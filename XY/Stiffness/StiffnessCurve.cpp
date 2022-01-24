#include "../TrigonalXYModel.cpp"
#include "../SquareXYModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

std::vector<double> twist_sampler(SquareXYModel *model) {
    return model->twist_stiffness();
}

int main(int argc, char* argv[]) {
    std::srand((unsigned)std::time( NULL ));

    std::string filename = argv[1];
    int N = std::stoi(argv[2]);
    std::cout << "N = " << N << std::endl;
    int num_threads = 4;

    const int L = 1;
    const float J = 1.;
    const float A = 0.;

    SquareXYModel *model = new SquareXYModel(N, L, J, 0., 0.);

    unsigned long long int resolution = 30;

#ifdef CLUSTER_UPDATE
    const int MCStep = 1;
#else
    const int MCStep = N*N*L;
#endif

    unsigned long long int steps_per_run = 10000*MCStep;

    unsigned long long int num_samples = 6000;
    unsigned long long int steps_per_sample = MCStep;


    std::vector<float> T(resolution);

    const float Tmax = 3.;
    const float Tmin = 0.1;
    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    auto start = std::chrono::high_resolution_clock::now();

    auto data = sample_r(twist_sampler, model, T, steps_per_run, num_samples, steps_per_sample, num_threads);
    string header = std::to_string(N) + "\t" + std::to_string(L);
    write_data(&data, &T, filename, header);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    std::cout << "Completion time: " << seconds/60. << " minutes." << std::endl;
}


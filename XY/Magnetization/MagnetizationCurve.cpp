#include "../SquareXYModel.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <chrono>

template <int q>
std::vector<float> magnetization_sampler(SquareXYModel *model) {
    return std::vector<float>{
                static_cast<float>(model->get_magnetization().norm()),
                static_cast<float>(model->energy()/model->V)
           };
}

int main(int argc, char* argv[]) {    
    int N = stoi(argv[2]);
    std::string filename = argv[1];
    int num_threads = 4;

    srand((unsigned)time( NULL ));

    cout << "N = " << N << endl;
    const int L = 1;
    const float J = 1.;
    const int q = 6;

    const int MCStep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, 0., 0.);

    const float Tmax = 3.;
    const float Tmin = 0.1;
    int resolution = 60;

    std::vector<float> T(resolution);

    for (int i = 0; i < resolution; i++) {
        T[i] = float(i)/resolution*Tmax + float(resolution - i)/resolution*Tmin;
    }

    int steps_per_run = 1000;//*MCStep;

    int num_samples = 500;
    int steps_per_sample = 1;//MCStep;

    unsigned long long nsteps = (unsigned long long) resolution*(steps_per_run + num_samples*steps_per_sample);

    auto start = std::chrono::high_resolution_clock::now();

    auto data = sample_r(magnetization_sampler<q>, model, &T, 10, steps_per_run, num_samples, steps_per_sample, num_threads);
    auto stats = summary_statistics(&data);
    std::string header = to_string(N) + "\t" + to_string(L);
    write_data(&stats, &T, filename, header); 

    auto stop = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int microseconds = duration.count();

    std::cout << "Duration: " << microseconds/1e6 << std::endl;
    std::cout << "Steps/s/thread: " << nsteps/float(microseconds)*1e6/num_threads << std::endl;
}

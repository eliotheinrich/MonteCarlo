#include "../SquareIsingModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

template <class ModelType>
void magnetization_run(ModelType *model, std::vector<float> *T, int num_runs,
                                                           int steps_per_run, 
                                                           int num_samples, 
                                                           int steps_per_sample,
                                                           int num_threads,
                                                           std::string filename) {

    int resolution = T->size();

    std::vector<MonteCarlo<ModelType>*> models(resolution);

    for (int i = 0; i < resolution; i++) {
        models[i] = new MonteCarlo<ModelType>(model->clone());
        models[i]->model->randomize_spins();
    }

    std::vector<std::vector<float>> M = std::vector<std::vector<float>>(resolution, std::vector<float>(num_samples*num_runs));
    std::vector<std::vector<float>> E = std::vector<std::vector<float>>(resolution, std::vector<float>(num_samples*num_runs));

    ctpl::thread_pool threads(num_threads);

    std::vector<std::future<void>> results(resolution);

    auto magnetization_samples = [num_runs, steps_per_run, num_samples, steps_per_sample](int id, int i, 
                                                                                          MonteCarlo<ModelType> *m, float T, 
                                                                                          std::vector<float> *M, 
                                                                                          std::vector<float> *E) {
        for (int n = 0; n < num_runs; n++) {
            m->model->randomize_spins();
            m->steps(steps_per_run, T);
            for (int i = 0; i < num_samples; i++) {
                (*M)[n*num_samples + i] = abs(m->model->get_magnetization());
                (*E)[n*num_samples + i] = m->model->energy()/m->model->V;
                m->steps(steps_per_sample, T);
            }
        }
    };

    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(magnetization_samples, i, models[i], (*T)[i], &M[i], &E[i]);
    }

    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }

    std::vector<float> avg_M(resolution);
    std::vector<float> avg_E(resolution);

    std::vector<float> err_M(resolution);
    std::vector<float> err_E(resolution);

    for (int i = 0; i < resolution; i++) {
        avg_M[i] = avg(&M[i]);
        avg_E[i] = avg(&E[i]);
        err_M[i] = stdev(&M[i], avg_M[i]);
        err_E[i] = stdev(&E[i], avg_E[i]);
    }

    std::ofstream output_file(filename);

    // Write header
    output_file << resolution << std::endl;

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t" << avg_M[i] << "\t" << err_M[i] << "\t" << avg_E[i] << "\t" << err_E[i] << std::endl;
    }

    output_file.close();
}


int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    int N = 32;
    int L = 1;
    const int MCStep = N*N*L;

    float J = 1.0;
    float B = 0;
    SquareIsingModel *model = new SquareIsingModel(N, L, J, B);

    int resolution = 50;
    int steps_per_run = 10000*MCStep;
    int num_samples = 100;
    int steps_per_sample = 50*MCStep;
    int num_runs = 1;
    int num_threads = 4;

    float T_max = 4.;
    float T_min = 0.05;

    std::vector<float> T(resolution);
    for (int i = 0; i < resolution; i++) {
        T[i] = T_max*i/resolution + T_min*(resolution - i)/resolution;
    }

    unsigned long long int nsteps = 3*resolution*(steps_per_run + num_samples*steps_per_sample);

    std::cout << "Number steps: " << nsteps << std::endl;
    std::cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << std::endl;

    auto start = std::chrono::high_resolution_clock::now();


    magnetization_run(model, &T, num_runs, steps_per_run, num_samples, steps_per_sample, num_threads, "MagnetizationCurve.txt");

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    std::cout << "Completion time: " << seconds/60. << " minutes." << std::endl;


}


#ifndef MONTECARLO
#define MONTECARLO

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctpl.h>
#include <thread>
#include <cmath>
#include "Utility.cpp"
#include "MonteCarlo.h"
#include "DataFrame.h"

#define PI 3.14159265
#define BOLTZMANN_CONSTANT 0.08617

double const_T(int n, int n_max, double Ti, double Tf) {
    return 0.5*(Tf + Ti);
}

double trig_T(int n, int n_max, double Ti, double Tf) {
    return Tf + 0.5*(Ti - Tf)*(1 - cos(n*PI/n_max));
}

double linear_T(int n, int n_max, double Ti, double Tf) {
    return Ti + (Tf - Ti)*(n_max - n)/double(n_max);
}

MonteCarlo::MonteCarlo(MCModel *model, int num_models) {
    this->num_models = 0;
    for (int i = 0; i < num_models; i++) {
        this->add_model(model->clone());
    }
}

MonteCarlo::MonteCarlo(std::vector<MCModel*> models) { 
    for (int i = 0; i < models.size(); i++) {
        this->add_model(models[i]);
    }
}

void MonteCarlo::add_model(MCModel *model) {
    this->models.push_back(model);
    this->random_generators.push_back(new std::minstd_rand()); // TODO seed?
    this->num_models++;
}

void MonteCarlo::steps(MCModel *model, unsigned long long nsteps, std::minstd_rand *rng) {
    // Performs MC simulation
    // nsteps: number of MC steps to perform

    double rf;
    double dE;
    double T = model->T;

    for (unsigned long long i = 0; i < nsteps; i++) {
        model->generate_mutation();
        dE = model->energy_change();

        rf = double((*rng)())/double(RAND_MAX);
        if (rf < std::exp(-dE/T)) {
            model->accept_mutation();
        } else {
            model->reject_mutation();
        }
    }
}

template<class dtype> 
DataFrame generate_samples(std::map<std::string, dtype> sampling_funcs(MCModel*), 
                            unsigned long long equilibration_steps,
                            unsigned long long num_samples,
                            unsigned long long steps_per_sample,
                            int num_threads) {
    ctpl::thread_pool threads(num_threads);

    std::vector<std::future<void>> results(num_models);

    // TODO annealing
    auto do_steps = [equilibration_steps](int id, MCModel *m) {
        MonteCarlo::steps(equilibration_steps);
    };

    // Do initial steps
    for (int i = 0; i < num_models; i++) {
        results[i] = threads.push(do_steps, models[i]);
    }

    // Join threads
    for (int i = 0; i < num_models; i++) {
        results[i].get();
    }

    std::vector<DataSlide> slides = std::vector<DataSlide>(0);


    auto take_samples = [num_samples, steps_per_sample, dtype_size, sampling_func](int id, MCModel *model, DataSlide *slide) {
        // Save params
        std::map<std::string, dtype> sample = sampling_func(model);
        // Initialize all keys into slide
        for (int n = 1; n < num_samples; n++) {
            sample = sampling_func(m->model);
            // Take samples
            MonteCarlo::steps(model, steps_per_sample);
        }

    };

    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(take_samples, i, models[i], &slides[i]);
    }

    // Join threads
    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }

    return arr;

}

#endif

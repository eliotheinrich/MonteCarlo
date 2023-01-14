#include "MonteCarlo.h"

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

std::map<std::string, Sample> to_sample(std::map<std::string, double> *s) {
    std::map<std::string, Sample> new_s;
    for (std::map<std::string, double>::iterator it = s->begin(); it != s->end(); ++it) {
        new_s.emplace(it->first, Sample(it->second));
    }
    return new_s;
};

DataFrame MonteCarlo::generate_samples(std::map<std::string, double> sampling_func(MCModel*), 
                                       unsigned long long equilibration_steps,
                                       unsigned long long num_samples,
                                       unsigned long long steps_per_sample,
                                       int num_threads,
                                       bool average_samples=true) {

    ctpl::thread_pool threads(num_threads);

    std::vector<std::future<void>> results(this->num_models);

    // TODO annealing
    auto do_steps = [equilibration_steps](int id, MCModel *model, std::minstd_rand *rng) {
        MonteCarlo::steps(model, equilibration_steps, rng);
    };

    // Do initial steps
    for (int i = 0; i < num_models; i++) {
        results[i] = threads.push(do_steps, models[i], random_generators[i]);
    }

    // Join threads
    for (int i = 0; i < num_models; i++) {
        results[i].get();
    }

    std::vector<DataSlide> slides = std::vector<DataSlide>(0);


    auto take_samples = [num_samples, steps_per_sample, sampling_func, average_samples](int id, MCModel *model, std::minstd_rand *rng, DataSlide *slide) {

        // Save params
        std::map<std::string, int> int_params = model->get_int_params();
        for (std::map<std::string, int>::iterator it = int_params.begin(); it != int_params.end(); ++it) {
            slide->add_int(it->first, it->second);
        }
        std::map<std::string, double> double_params = model->get_double_params();
        for (std::map<std::string, double>::iterator it = double_params.begin(); it != double_params.end(); ++it) {
            slide->add_double(it->first, it->second);
        }

        // Initialize all keys from sampling funcs
        std::map<std::string, double> raw_samples = sampling_func(model);
        for (std::map<std::string, double>::iterator it = raw_samples.begin(); it != raw_samples.end(); ++it) {
            slide->add_data(it->first); 
        }

        std::map<std::string, Sample> samples = to_sample(&raw_samples);

        for (int n = 0; n < num_samples; n++) {
            // Take samples
            MonteCarlo::steps(model, steps_per_sample, rng);

            raw_samples = sampling_func(model);
            std::map<std::string, Sample> new_samples = to_sample(&raw_samples);
            if (average_samples) {
                for (std::map<std::string, Sample>::iterator it = new_samples.begin(); it != new_samples.end(); ++it) {
                    samples[it->first] = samples[it->first].combine(&new_samples[it->first]);
                }
            } else {
                for (std::map<std::string, Sample>::iterator it = new_samples.begin(); it != new_samples.end(); ++it) {
                    slide->push_data(it->first, it->second);
                }
            }
        }
    };

    for (int i = 0; i < num_models; i++) {
        results[i] = threads.push(take_samples, models[i], random_generators[i], &slides[i]);
    }

    // Join threads
    for (int i = 0; i < num_models; i++) {
        results[i].get();
    }

    return DataFrame(slides);

}
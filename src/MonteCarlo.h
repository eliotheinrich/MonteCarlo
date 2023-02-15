#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "ctpl.h"
#include <thread>
#include <cmath>
#include <map>
#include "Utility.h"
#include "DataFrame.h"

#define PI 3.14159265
#define BOLTZMANN_CONSTANT 0.08617

enum cooling_type {
    constant,
    cosine,
    linear,
    exponential
};

inline double const_T(int n, int n_max, double Ti, double Tf) {
    return 0.5*(Tf + Ti);
}

inline double trig_T(int n, int n_max, double Ti, double Tf) {
    return Tf + 0.5*(Ti - Tf)*(1 - cos(n*PI/n_max));
}

inline double linear_T(int n, int n_max, double Ti, double Tf) {
    return Ti + (Tf - Ti)*(n_max - n)/double(n_max);
}

class MCModel {
    // Most basic Monte-Carlo model to be simulated must have some notion of energy
    // as well as a mutation data structure. Specifics must be supplied by child classes.
    public:
        double T;
        double stored_energy;
        unsigned long long nsteps;

        virtual double energy() const = 0;
        virtual double energy_change() = 0;
        virtual void generate_mutation() = 0;
        virtual void accept_mutation() = 0;
        virtual void reject_mutation() = 0;
        virtual MCModel* clone() = 0;
        virtual std::map<std::string, double> get_double_params() const {
            return std::map<std::string, double>();
        }
        virtual std::map<std::string, int> get_int_params() const {
            return std::map<std::string, int>();
        }
};

template<class ModelType>
class MonteCarlo {
    // MonteCarlo is a wrapper for a MCModel
    // Equipped with tools for doing Monte-Carlo simulation on generic systems
    // with well-defined energy.

    private:
        void init();
        std::vector<ModelType*> models;
        int num_models;

        // Random number generators
        std::vector<std::minstd_rand*> random_generators;

    public: 
        MonteCarlo(ModelType *model, int num_models) {
            this->num_models = 0;
            for (int i = 0; i < num_models; i++) {
                this->add_model(model->clone());
            }
        }

        MonteCarlo(std::vector<ModelType*> models) {
            this->num_models = 0;
            for (int i = 0; i < models.size(); i++) {
                this->add_model(models[i]);
            }
        }

        inline int get_num_models() {
            return num_models;
        }

        void add_model(ModelType *model) {
            this->models.push_back(model);
            this->random_generators.push_back(new std::minstd_rand()); // TODO seed?
            this->num_models++;
        }

        static void steps(ModelType *model, unsigned long long num_steps, std::minstd_rand* r) {
            // Performs MC simulation
            // nsteps: number of MC steps to perform

            double rf;
            double dE;
            double T = model->T;

            for (unsigned long long i = 0; i < num_steps; i++) {
                model->generate_mutation();
                dE = model->energy_change();

                rf = double((*r)())/double(RAND_MAX);
                if (rf < std::exp(-dE/T)) {
                    model->accept_mutation();
                } else {
                    model->reject_mutation();
                }
            }
        }

        DataFrame generate_samples(std::map<std::string, double> sampling_func(ModelType*), 
                                   unsigned long long equilibration_steps,
                                   unsigned long long num_samples,
                                   unsigned long long steps_per_sample,
                                   int num_threads,
                                   bool average_samples=true) {
            ctpl::thread_pool threads(num_threads);

            std::vector<std::future<void>> results(this->num_models);

            // TODO annealing
            auto do_steps = [equilibration_steps](int id, ModelType *model, std::minstd_rand *rng) {
                MonteCarlo<ModelType>::steps(model, equilibration_steps, rng);
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
            for (int i = 0; i < num_models; i++) {
                slides.push_back(DataSlide());
            }

            auto take_samples = [num_samples, steps_per_sample, sampling_func, average_samples](int id, ModelType *model, std::minstd_rand *rng, DataSlide *slide) {
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
                std::vector<std::string> keys;
                for (std::map<std::string, double>::iterator it = raw_samples.begin(); it != raw_samples.end(); ++it) {
                    slide->add_data(it->first); 
                    keys.push_back(it->first);
                }

                std::map<std::string, Sample> samples = to_sample(&raw_samples);

                for (int n = 0; n < num_samples-1; n++) {
                    // Take samples
                    MonteCarlo<ModelType>::steps(model, steps_per_sample, rng);

                    raw_samples = sampling_func(model);
                    std::map<std::string, Sample> new_samples = to_sample(&raw_samples);
                    if (average_samples) {
                        for (const auto &key : keys) {
                            samples[key] = samples[key].combine(&new_samples[key]);
                        }
                    } else {
                        for (std::map<std::string, Sample>::iterator it = new_samples.begin(); it != new_samples.end(); ++it) {
                            slide->push_data(it->first, it->second);
                        }
                    }
                }

                if (average_samples) {
                    for (std::map<std::string, Sample>::iterator it = samples.begin(); it != samples.end(); ++it) {
                        slide->push_data(it->first, it->second);
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
};

#endif
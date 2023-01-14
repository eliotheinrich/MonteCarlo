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
#include "Utility.cpp"
#include "DataFrame.h"

#define PI 3.14159265
#define BOLTZMANN_CONSTANT 0.08617

enum cooling_type {
    constant,
    cosine,
    linear,
    exponential
};

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
        virtual std::map<std::string, double> get_double_params() const;
        virtual std::map<std::string, int> get_int_params() const;
};

class MonteCarlo {
    // MonteCarlo is a wrapper for a MCModel
    // Equipped with tools for doing Monte-Carlo simulation on generic systems
    // with well-defined energy.

    private:
        void init();
        std::vector<MCModel*> models;
        int num_models;

        // Random number generators
        std::vector<std::minstd_rand*> random_generators;

    public: 
        MonteCarlo(MCModel *model, int num_models);
        MonteCarlo(std::vector<MCModel*> models);

        void add_model(MCModel *model);

        static void steps(MCModel *model, unsigned long long num_steps, std::minstd_rand* r);

        DataFrame generate_samples(std::map<std::string, double> sampling_func(MCModel*), 
                                   unsigned long long equilibration_steps,
                                   unsigned long long num_samples,
                                   unsigned long long steps_per_sample,
                                   int num_threads,
                                   bool average_samples);
};

#endif
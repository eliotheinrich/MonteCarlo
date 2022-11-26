#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctpl.h>
#include <thread>
#include <cmath>
#include <map>
#include "Utility.cpp"

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
        virtual double energy_change() const = 0;
        virtual void generate_mutation() = 0;
        virtual void accept_mutation() = 0;
        virtual void reject_mutation() = 0;
        virtual MCModel* clone() = 0;
        virtual std::map<std::string, double> get_double_params() const;
        virtual std::map<std::string, int> get_int_params() const;
};

class DataSampler {
    public:
        virtual std::vector<double> collect_data(const MCModel*) { return std::vector<double>(0); }
};

class DataFrame {
    public:
        std::map<std::string, std::vector<std::pair<double, double>>> data;

        std::vector<std::pair<double, double>>& operator[](std::string s);
        void emplace(std::string s, std::vector<std::pair<double, double>> p);
        int size();
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

        std::vector<DataFrame*> data;

    public:
        
        MonteCarlo(MCModel *model, int num_models);
        MonteCarlo(std::vector<MCModel*> models);

        void step_all_models(unsigned long long num_steps) {

        } 
        static void steps(MCModel* model, unsigned long long num_steps, std::minstd_rand* r);
        static void cooling_steps(MCModel* model, 
                                  unsigned long long num_steps, 
                                  int num_updates, 
                                  double Ti, 
                                  double Tf,
                                  std::minstd_rand* r, 
                                  cooling_type cs);

        static void collect_data(MCModel* model, std::map<std::string, DataSampler*> samplers, DataFrame* data);
        void sample(std::map<std::string, DataSampler*> samplers, 
                    int num_cooling_steps,
                    int steps_per_sample, 
                    int num_samples,
                    int num_threads,
                    double Ti,
                    double Tf);

        bool write_data(std::string filename);
};

#endif
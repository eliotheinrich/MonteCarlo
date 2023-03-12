#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <cmath>
#include <map>
#include <random>
#include <string>
#include <DataFrame.hpp>
#include <Simulator.hpp>

#define PI 3.14159265
//#define BOLTZMANN_CONSTANT 0.08617

#define DEFAULT_COOLING_SCHEDULE "constant"
#define DEFAULT_NUM_COOLING_STEPS 100
#define DEFAULT_RANDOM_SEED -1

enum CoolingSchedule {
    Constant,
    Trig,
    Linear,
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

static CoolingSchedule parse_cooling_schedule(std::string s) {
    if ((s == "constant") || (s == "const")) { return CoolingSchedule::Constant; }
    else if ((s == "trig") || (s == "cosine")) { return CoolingSchedule::Trig; }
    else if ((s == "linear")) { return CoolingSchedule::Linear; }
    else { return CoolingSchedule::Constant; }
}

class MCModel {
    // Most basic Monte-Carlo model to be simulated must have some notion of energy
    // as well as a mutation data structure. Specifics must be supplied by child classes.
    private:
        std::minstd_rand *rng;

    public:
        double temperature;


        // Give all MCModels a source of random numbers
        int rand() { return (*rng)(); }
        float randf() { return float(rand())/float(RAND_MAX); }
        MCModel(int seed=DEFAULT_RANDOM_SEED) : rng(new std::minstd_rand(seed)) {}

        // To be overridden by child classes
        virtual double energy() const = 0;
        virtual double energy_change() = 0;
        virtual void generate_mutation() = 0;
        virtual void accept_mutation() = 0;
        virtual void reject_mutation() = 0;

        // Can optionally delay initialization dynamic variables until after processes have been forked by moving
        // expensive operations to the init function
        // This may lead to performance increases for small timesteps and large system sizes
        virtual void init() {}
        virtual MCModel* clone(Params &params) = 0;

        virtual std::map<std::string, Sample> take_samples() const {
            return std::map<std::string, Sample>();
        }
};

class MonteCarloSimulator : public Simulator {
    private:
        std::minstd_rand *rng;
        int random_seed;

        // For annealing
        CoolingSchedule cooling_schedule;
        uint num_cooling_steps;
        double init_temperature;
        double temperature;

        MCModel *model;
    
    public:
        MonteCarloSimulator(Params &params, MCModel *model);

        uint rand() { return (*rng)(); }
        float randf() { return float(rand())/float(RAND_MAX); }

        virtual void init_state();
        virtual void timesteps(uint num_steps);
        virtual void equilibration_timesteps(uint num_steps);
        virtual Simulator* clone(Params &params) { return new MonteCarloSimulator(params, model); }
        virtual std::map<std::string, Sample> take_samples() { return model->take_samples(); };
};

// TODO remove
static std::vector<std::string> split(std::string *s, std::string delim) {
    std::vector<std::string> vals(0);

    std::string str = *s;
    int pos = 0;
    std::string token;
    while ((pos = str.find(delim)) != std::string::npos) {
        token = str.substr(0, pos);
        vals.push_back(token);
        str.erase(0, pos + delim.length());
    }

    vals.push_back(str);

    return vals;
}

static const inline int mod(int a, int b) {
    int c = a % b;
    return (c < 0) ? c + b : c;
}

#endif
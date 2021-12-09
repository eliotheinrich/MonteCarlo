#ifndef MONTECARLO
#define MONTECARLO

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <cmath>
#include "Utility.cpp"

#define PI 3.14159265

using namespace std;

template <class MCModel>
class MonteCarlo;

class MCModel {
    // Most basic Monte-Carlo model to be simulated must have some notion of energy
    // as well as a mutation data structure. Specifics must be supplied by child classes.
    public:
        // Need to implement:
        virtual const float energy()=0;
        virtual const float energy_change()=0;
        virtual void generate_mutation()=0;
        virtual void accept_mutation()=0;
        virtual void reject_mutation()=0;
        virtual MCModel* clone()=0;
};

float const_T(int n, int n_max, float T_min, float T_max) {
    return 0.5*(T_max + T_min);
}

float trig_T(int n, int n_max, float T_min, float T_max) {
    return T_max + 0.5*(T_min - T_max)*(1 - cos(n*PI/n_max));
}

float linear_T(int n, int n_max, float T_min, float T_max) {
    return T_min + (T_max - T_min)*(n_max - n)/float(n_max);
}

template<class MCModel>
class MonteCarlo {
    // MonteCarlo is a wrapper for a MCModel
    // Equipped with tools for doing Monte-Carlo simulation on generic systems
    // with well-defined energy.

    public:
        MCModel *model;
        int accepted;
        long nsteps;
        float energy;

        MonteCarlo(MCModel *model) {
            this->model = model;
            this->accepted = 0.;
            this->nsteps = 0;
            this->energy = model->energy();
        }

        template<class A>
        vector<A> sample(function<A(MCModel*)> f, float T, int num_samples, int steps_per_sample) {
            vector<A> samples(num_samples);
            for (int i = 0; i < num_samples; i++) {
                samples[i] = f(model);
                steps(steps_per_sample, T);
            }

            return samples;
        }


        void steps(int nsteps, float T) {
            // Performs MC simulation
            // nsteps: number of MC steps to perform
            // T: temperature

            float r;
            float dE;

            for (int i = 0; i < nsteps; i++) {
                model->generate_mutation();
                dE = model->energy_change();

                r = float(rand())/float(RAND_MAX);
                if (r < exp(-dE/T)) {
                    accepted++;

                    model->accept_mutation();
                    this->energy += dE;
                } else {
                    model->reject_mutation();
                }
            }


            this->nsteps += nsteps;
        }
};

// run_MC with everything
template <class MCModel, class A>
vector<A> run_MC(MonteCarlo<MCModel> *m, int nsteps, string cooling,
                                                     float T_max,
                                                     float T_min,
                                                     int num_updates, 
                                                     function<A(MCModel*)> f) {

    // Establish cooling schedule
    float (*update_T)(int n, int n_max, float T_min, float T_max);
    if (cooling == "auto") {
        T_min = T_max;
        update_T = *const_T;
    } else {
        if (T_max == -1 || T_min == -1) {
            cout << "Need to supply T_max and T_min!" << endl;
        }
        if (cooling == "trig") {
            update_T = *trig_T;
        } else if (cooling == "linear") {
            update_T = *linear_T;
        } else if (cooling == "const") {
            update_T = *const_T;
        }
    }
    
    vector<A> log(num_updates);
    int update_freq = nsteps/num_updates;

    float T;
    for (int i = 0; i < num_updates; i++) {
        T = update_T(i, num_updates, T_min, T_max);
        m->steps(update_freq, T); 

        log[i] = f(m->model);
    }

    return log;
}


// Basic run_MC
template <class MCModel>
void run_MC(MonteCarlo<MCModel> *m, int nsteps, float T) {
    m->steps(nsteps, T);
}

// run_MC with cooling schedule
template <class MCModel>
void run_MC(MonteCarlo<MCModel> *m, int nsteps, string cooling, float Tmin, float Tmax, int num_updates = 100) {
    run_MC(m, nsteps, cooling, Tmin, Tmax, num_updates, function<int(MCModel*)>([](MCModel* g) { return 0; }));
}

// run_MC with logging function
template <class MCModel, class A>
vector<A> run_MC(MonteCarlo<MCModel> *m, int nsteps, float T, function<A(MCModel*)> f, int num_logitems) {
    return run_MC(m, nsteps, "const", T, T, num_logitems, f);
}

template <class MCModel>
vector<MonteCarlo<MCModel>*> parallel_tempering(MCModel *model, vector<float> Ts, int steps_per_exchange, 
                                                                                  int num_exchanges, 
                                                                                  int equilibration_steps = -1) {

    int num_threads = Ts.size();
    vector<thread> threads(num_threads);
    vector<MonteCarlo<MCModel>*> models(num_threads);

    // Initialize models
    for (int i = 0; i < num_threads; i++) {
        models[i] = new MonteCarlo<MCModel>(model->clone());
    }

    int accepted = 0;
    int rejected = 0;

    float r;
    float dE;
    float dB;
    MonteCarlo<MCModel> *model_buffer;

    for (int k = 0; k < num_exchanges; k++) {
        // Give threads work
        for (int i = 0; i < num_threads; i++) {
            threads[i] = thread(&MonteCarlo<MCModel>::steps, models[i], steps_per_exchange, Ts[i]);
        }

        // Join threads
        for (int i = 0; i < num_threads; i++) {
            threads[i].join();
        }

        // Make exchanges
        for (int i = 0; i < num_threads-1; i++) {
            r = float(rand())/float(RAND_MAX);
            dE = models[i]->energy - models[i+1]->energy;
            dB = 1./Ts[i] - 1./Ts[i+1];
            if (r < exp(dE*dB)) {
                accepted++;
                model_buffer = models[i];
                models[i] = models[i+1];
                models[i+1] = model_buffer;
            } else {
                rejected++;
            }
        }
    }

    // Final equilibration
    if (equilibration_steps != -1) {
        for (int i = 0; i < num_threads; i++) {
            threads[i] = thread(&MonteCarlo<MCModel>::steps, models[i], equilibration_steps, Ts[i]);
        }

        for (int i = 0; i < num_threads; i++) {
            threads[i].join();
        }
    }

    return models;
}


#endif

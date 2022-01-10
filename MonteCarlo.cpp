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

#define PI 3.14159265
#define BOLTZMANN_CONSTANT 0.08617


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

float const_T(int n, int n_max, float Ti, float Tf) {
    return 0.5*(Tf + Ti);
}

float trig_T(int n, int n_max, float Ti, float Tf) {
    return Tf + 0.5*(Ti - Tf)*(1 - cos(n*PI/n_max));
}

float linear_T(int n, int n_max, float Ti, float Tf) {
    return Ti + (Tf - Ti)*(n_max - n)/float(n_max);
}

template<class MCModel>
class MonteCarlo {
    // MonteCarlo is a wrapper for a MCModel
    // Equipped with tools for doing Monte-Carlo simulation on generic systems
    // with well-defined energy.

    public:
        MCModel *model;
        minstd_rand r;
        int accepted;
        unsigned long long nsteps;
        float energy;

        MonteCarlo(MCModel *model) {
            this->model = model;
            this->accepted = 0.;
            this->nsteps = 0;
            this->energy = model->energy();
            this->r.seed(rand());
        }

        void steps(unsigned long long nsteps, float T) {
            // Performs MC simulation
            // nsteps: number of MC steps to perform
            // T: temperature

            float rf;
            float dE;

            for (unsigned long long i = 0; i < nsteps; i++) {
                model->generate_mutation();
                dE = model->energy_change();
                
                rf = float(r())/float(RAND_MAX);
                if (rf < exp(-dE/T)) {
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
vector<A> run_MC(MonteCarlo<MCModel> *m, unsigned long long nsteps, string cooling,
                                                                    float Ti,
                                                                    float Tf,
                                                                    int num_updates, 
                                                                    function<A(MCModel*)> f) {

    // Establish cooling schedule
    float (*update_T)(int n, int n_max, float Ti, float Tf);
    if (cooling == "auto") {
        Ti = Tf;
        update_T = *const_T;
    } else {
        if (Ti == -1 || Tf == -1) {
            cout << "Need to supply Ti and Tf!" << endl;
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
        T = update_T(i, num_updates, Ti, Tf);
        m->steps(update_freq, T); 

        log[i] = f(m->model);
    }

    return log;
}


// Basic run_MC
template <class MCModel>
void run_MC(MonteCarlo<MCModel> *m, unsigned long long nsteps, float T) {
    m->steps(nsteps, T);
}

// run_MC with cooling schedule
template <class MCModel>
void run_MC(MonteCarlo<MCModel> *m, unsigned long long nsteps, string cooling, float Ti, float Tf, int num_updates = 100) {
    run_MC(m, nsteps, cooling, Ti, Tf, num_updates, function<int(MCModel*)>([](MCModel* g) { return 0; }));
}

// run_MC with logging function
template <class MCModel, class A>
vector<A> run_MC(MonteCarlo<MCModel> *m, unsigned long long nsteps, float T, function<A(MCModel*)> f, int num_logitems) {
    return run_MC(m, nsteps, "const", T, T, num_logitems, f);
}


#endif

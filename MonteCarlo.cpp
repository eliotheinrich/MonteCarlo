#ifndef MONTECARLO
#define MONTECARLO

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

template <class MCModel>
class LogItem {
    // Stores some datatype specified by child classes
    public:
        LogItem() {}
        LogItem(MonteCarlo<MCModel> *m, MCModel *model) {}
        friend ostream& operator<<(ostream& os, const LogItem& logitem);
};

template <class LogItemType>
class Log {
    // Stores an array of some LogItem updated at some specified frequency in simulation
    public:
        int len;
        vector<LogItemType> items;

        Log() {
            this->len = 0;
            items = vector<LogItemType>(0);
        }

        Log(int len) {
            this->len = len;
            items = vector<LogItemType>(len);
        }

        LogItemType & operator[](int i) {
            return items[i];
        }

        void save_log(string filename) {
            ofstream output_file;
            output_file.open(filename);
            for (int i = 0; i<len; i++) {
                    output_file << items[i] << endl; 
            }
            output_file.close();
        }
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
        int step;
        float energy;

        MonteCarlo(MCModel *model) {
            this->model = model;
            this->accepted = 0.;
            this->step = 0;
            this->energy = model->energy();
        }

        template<class A>
        A expectation(function<A(MCModel*)> f, float T, int num_samples, int steps_per_sample) {
            A avg_A = 0.;
            for (int i = 0; i < num_samples; i++) {
                avg_A += f(model)/num_samples;
                steps(steps_per_sample, T);
            }

            return avg_A;
        }


        void steps(int nsteps, float T) {
            // Performs MC simulation
            // nsteps: number of MC steps to perform
            // T: temperature

            float r;
            float dE;

            for (int step = 0; step<nsteps; step++) {

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

            this->step += nsteps;
        }
};

// run_MC with everything
template <class MCModel, class A>
void run_MC(MonteCarlo<MCModel> *m, int nsteps, string cooling,
                                                float T_max,
                                                float T_min,
                                                int num_updates, 
                                                bool record_log,
                                                string filename,
                                                int num_logitems, 
                                                function<A(MCModel*)> f) {

    // Establish cooling schedule
    float (*update_T)(int n, int n_max, float T_min, float T_max);
    if (cooling == "auto") {
        T_min = T_max;
        update_T = *const_T;
    } else {
        if (T_max == -1 || T_min == -1) {
            cout << "Need to supply T_max and T_min!" << endl;
            return;
        }
        if (cooling == "trig") {
            update_T = *trig_T;
        } else if (cooling == "linear") {
            update_T = *linear_T;
        } else if (cooling == "const") {
            update_T = *const_T;
        }
    }

    int update_freq = nsteps/num_updates;

    float T;
    for (int i = 0; i < num_updates; i++) {
        T = update_T(i, num_updates, T_min, T_max);
        m->steps(update_freq, T); 

    }
}

// Basic run_MC
template <class MCModel>
void run_MC(MonteCarlo<MCModel> *m, int nsteps, float T) {
    m->steps(nsteps, T);
}

// run_MC with cooling schedule
template <class MCModel>
void run_MC(MonteCarlo<MCModel> *m, int nsteps, string cooling, float Tmin, float Tmax, int num_updates = 100) {
    run_MC(m, nsteps, cooling, Tmin, Tmax, num_updates, false, "", 0, function<void(MCModel*)>([](MCModel* g) {}));
}

// run_MC with logging function
template <class MCModel, class A>
void run_MC(MonteCarlo<MCModel> *m, int nsteps, float T, int num_logitems, function<A(MCModel*)> f, string filename = "Log.txt") {
    run_MC(m, nsteps, "const", T, T, 1, true, filename, num_logitems, f);
}

template <class MCModel>
vector<MCModel*> parallel_tempering(MCModel *model, vector<float> Ts, int steps_per_exchange, int num_exchanges) {
    int num_threads = Ts.size();
    vector<thread> threads(num_threads);
    vector<MCModel*> models(num_threads);
    vector<MonteCarlo<MCModel>*> mc_models(num_threads);

    // Initialize models
    for (int i = 0; i < num_threads; i++) {
        models[i] = model->clone();
        mc_models[i] = new MonteCarlo<MCModel>(models[i]);
    }

    int accepted = 0;
    int rejected = 0;

    float r;
    float dE;
    float dB;
    MonteCarlo<MCModel> *mc_model_buffer;
    MCModel *model_buffer;
    for (int k = 0; k < num_exchanges; k++) {
        // Give threads work
        for (int i = 0; i < num_threads; i++) {
            threads[i] = thread(&MonteCarlo<MCModel>::steps, mc_models[i], steps_per_exchange, Ts[i]);
        }

        // Join threads
        for (int i = 0; i < num_threads; i++) {
            threads[i].join();
        }

        // Make exchanges
        for (int i = 0; i < num_threads-1; i++) {
            r = float(rand())/float(RAND_MAX);
            dE = mc_models[i]->energy - mc_models[i+1]->energy;
            dB = 1./Ts[i] - 1./Ts[i+1];
            //cout << "dE = " << dE << ", dB = " << dB << ", dB*dE = " << dB*dE << endl;
            if (r < exp(dE*dB)) {
                //cout << "Made a swap between sites " << i << " and " << i+1 << endl;
                accepted++;
                mc_model_buffer = mc_models[i];
                mc_models[i] = mc_models[i+1];
                mc_models[i+1] = mc_model_buffer;

                model_buffer = models[i];
                models[i] = models[i+1];
                models[i+1] = model_buffer;
            } else {
                rejected++;
                //cout << "Swap rejected between " << i << " and " << i+1 << endl;
            }
        }
    }

    //cout << "Accepted: " << accepted << endl;
    //cout << "Rejected: " << rejected << endl;
    //cout << "Ratio: " << float(accepted)/(rejected + accepted) << endl;

    // Final equilibration
    for (int i = 0; i < num_threads; i++) {
        cout << mc_models[i]->energy << endl;
        cout << models[i]->energy() << endl;
    }
    cout << endl;

    for (int i = 0; i < num_threads; i++) {
        threads[i] = thread(&MonteCarlo<MCModel>::steps, mc_models[i], 5*steps_per_exchange, Ts[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        threads[i].join();
        cout << models[i]->energy() << endl;
    }
    
    return models;
}


#endif

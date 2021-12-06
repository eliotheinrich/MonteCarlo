#ifndef MONTECARLO
#define MONTECARLO

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "Utility.cpp"

#define PI 3.14159265

using namespace std;

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
};

class LogItem {
    // Stores some datatype specified by child classes
    public:
        LogItem() {}
        LogItem(MonteCarlo *m, MCModel *model) {}
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

class EnergyLogItem : public LogItem {
    public:
        float energy;

        EnergyLogItem() {}

        EnergyLogItem(MonteCarlo *m, MCModel *model) {
            this->energy = m->energy;
        }

        friend ostream& operator<<(ostream& os, const EnergyLogItem& logitem) {
            os << logitem.energy;
            return os;
        }
};

//template <typename T> string type_name();

template <class LogItemType=EnergyLogItem, class Model>
void run_MC(Model *model, int nsteps, string cooling = "auto",
                                      float T_max = -1,
                                      float T_min = -1,
                                      bool record_log = false,
                                      int num_updates = 100, 
                                      string filename = "Log.txt") {

    MonteCarlo *m = new MonteCarlo(model);

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

    Log<LogItemType> log = Log<LogItemType>(num_updates);

    int update_freq = nsteps/num_updates;

    float T;
    for (int i = 0; i < num_updates; i++) {
        T = update_T(i, num_updates, T_min, T_max);
        m->steps(update_freq, T); 

        if (record_log) {
            log[i] = LogItemType(m, model);
        }
    }

    if (record_log) { log.save_log(filename); }
}

#endif

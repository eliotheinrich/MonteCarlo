#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include <iostream>
#include <math.h>

#define BOLTZMANN_CONSTANT 0.08617

using namespace std;
using namespace Eigen;

vector<vector<float>> take_data(TrigonalModel *model, vector<float> *Ts, int equilibration_steps, 
                                                                         int num_samples, 
                                                                         int steps_per_sample,
                                                                         int num_threads) {

    int resolution = Ts->size();
    model->randomize_spins();
    MonteCarlo<TrigonalModel> *m = new MonteCarlo<TrigonalModel>(model);

    vector<MonteCarlo<TrigonalModel>*> ms(resolution);

    for (int i = 0; i < resolution; i++) {
        ms[i] = new MonteCarlo<TrigonalModel>(model->clone());
    }

    vector<vector<float>> X = vector<vector<float>>(resolution, vector<float>(num_samples));

    auto start = chrono::high_resolution_clock::now();

    auto susceptibility_samples = [equilibration_steps, num_samples, steps_per_sample](int id, int i, float T, 
                                                                                      MonteCarlo<TrigonalModel> *m, 
                                                                                      vector<float> *samples) {
        m->steps(equilibration_steps, T);
        for (int j = 0; j < num_samples; j++) {
            (*samples)[j] = m->model->get_magnetization().dot(m->model->B)/pow(m->model->B.norm(), 2);
            m->steps(steps_per_sample, T);
        }
    };

    ctpl::thread_pool threads(num_threads);
    vector<future<void>> results(resolution);
    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(susceptibility_samples, i, (*Ts)[i], ms[i], &X[i]);
    }

    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }
    

    return X;
}

void write_data(vector<vector<float>> X, vector<float> Ts, string filename) {
    ofstream output_file(filename);

    int resolution = X.size();
    int num_samples = X[0].size();

    // Write header
    output_file << resolution << "\t" << num_samples << endl;

    float avg_X;
    for (int i = 0; i < resolution; i++) {
        output_file << Ts[i]/BOLTZMANN_CONSTANT << "\t";
        for (int j = 0; j < num_samples; j++) {
            output_file << "(" << X[i][j] << ")";
            if (j < num_samples - 1) { output_file << '\t'; }
        }
        output_file << endl;
    }

    output_file.close();
}


int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    // Folder to store data
    string foldername = argv[1];

    // Load params config file
    ifstream paramfile(argv[2]);

    float S = 3.5;

    // Magnetic field strength
    float Bm = S*stof(argv[3]);

    // Load paramfile line by line
    string line;
    getline(paramfile, line);
    getline(paramfile, line);
    vector<string> paramss = split(&line, ",");

    int N = stoi(paramss[0]);
    int L = stoi(paramss[1]);
    const int MCStep = N*N*L;

    float J1 = S*S*stof(paramss[2]);
    float J2 = S*S*stof(paramss[3]);
    float K1 = S*stof(paramss[4]);
    float K2 = S*stof(paramss[5]);
    float K3 = S*stof(paramss[6]);

    getline(paramfile, line);
    getline(paramfile, line);
    getline(paramfile, line);
    paramss = split(&line, ",");

    int resolution = stoi(paramss[0]);
    int steps_per_run = stoi(paramss[1])*MCStep;
    int num_samples = stoi(paramss[2]);
    int steps_per_sample = stoi(paramss[3])*MCStep;
    int num_runs = stoi(paramss[4]);

    float T_max = 60*BOLTZMANN_CONSTANT; // In Kelvin
    float T_min = 0.1*BOLTZMANN_CONSTANT;
    vector<float> Ts(resolution);
    for (int i = 0; i < resolution; i++) {
        Ts[i] = T_max*i/resolution + T_min*(resolution - i)/resolution;
    }

    Vector3f Bhat1; Bhat1 << 0., 1., 0.;
    Vector3f Bhat2; Bhat2 << 1., 0., 0.;
    Vector3f Bhat3; Bhat3 << 0., 0., 1.;
    Vector3f B1 = Bm*Bhat1;
    Vector3f B2 = Bm*Bhat2; 
    Vector3f B3 = Bm*Bhat3;

    TrigonalModel *model1 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B1);
    TrigonalModel *model2 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B2);
    TrigonalModel *model3 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B3);

    vector<vector<float>> X1 = vector<vector<float>>(resolution, vector<float>(num_samples));
    vector<vector<float>> X2 = vector<vector<float>>(resolution, vector<float>(num_samples));
    vector<vector<float>> X3 = vector<vector<float>>(resolution, vector<float>(num_samples));

    int num_threads = 30;
    unsigned long long int nsteps = 3*resolution*num_runs*(steps_per_run + num_samples*steps_per_sample);

    cout << "Number steps: " << nsteps << endl;
    cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << endl;

    auto start = chrono::high_resolution_clock::now();


    for (int i = 0; i < num_runs; i++) {
        X1 = take_data(model1, &Ts, steps_per_run, num_samples, steps_per_sample, num_threads);
        X2 = take_data(model2, &Ts, steps_per_run, num_samples, steps_per_sample, num_threads);
        X3 = take_data(model3, &Ts, steps_per_run, num_samples, steps_per_sample, num_threads);

        write_data(X1, Ts, foldername + "/SusceptibilityCurve1_run" + to_string(i) + ".txt");
        write_data(X2, Ts, foldername + "/SusceptibilityCurve2_run" + to_string(i) + ".txt");
        write_data(X3, Ts, foldername + "/SusceptibilityCurve3_run" + to_string(i) + ".txt");
    }


    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    cout << "Completion time: " << seconds/60. << " minutes." << endl;


}


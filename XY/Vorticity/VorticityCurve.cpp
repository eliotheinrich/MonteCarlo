#include "../TrigonalXYModel.cpp"
#include "../SquareXYModel.cpp"
#include <ctpl_stl.h>
#include <iostream>
#include <chrono>

using namespace std;

void func(string filename, int N, int num_threads) {
    srand((unsigned)time( NULL ));

    cout << "N = " << N << endl;
    const int L = 1;
    const float J = 1.;
    const float B = 0.;
    const float Bp = 0.;

    const int MCStep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);
    //TrigonalXYModel *model = new TrigonalXYModel(N, L, J);
    MonteCarlo<SquareXYModel> *m = new MonteCarlo<SquareXYModel>(model);


    const float Tmax = 4.0;
    const float Tmin = 0.3;
    int res = 50;

    vector<float> Ts(res);
    ofstream output(filename);

    for (int i = 0; i < res; i++) {
        Ts[i] = float(i)/res*Tmax + float(res - i)/res*Tmin;
    }

    //int num_exchanges = 5000;
    //int steps_per_exchange = 5*MCStep;
    //int equilibration_steps = 5000*MCStep;

    int num_exchanges = 1000;
    int steps_per_exchange = 4*MCStep;
    int equilibration_steps = 1000*MCStep;
    auto models = parallel_tempering(model, Ts, num_exchanges, steps_per_exchange, equilibration_steps, num_threads);

    int num_samples = 3000;
    int steps_per_sample = 5*MCStep;
    ctpl::thread_pool threads(num_threads);
    vector<future<vector<vector<double>>>> results(res);
    vector<vector<double>> samples(num_samples);

    auto twist_sampling = [&models, &Ts, steps_per_sample, num_samples](int id, int i) {
        vector<vector<double>> samples(num_samples);
        for (int j = 0; j < num_samples; j++) {
            samples[j] = models[i]->model->vorticity();
            models[i]->steps(steps_per_sample, Ts[i]);
        }
        return samples;
    };

    // Write header
    output << res << "\t" << num_samples << endl;

    // Give threads jobs
    for (int i = 0; i < res; i++) {
        results[i] = threads.push(twist_sampling, i);
    } 

    for (int i = 0; i < res; i++) {
        samples = results[i].get();
        output << Ts[i] << "\t";
        for (int j = 0; j < num_samples; j++) {
            output << "(" << samples[j][0] << ", " << samples[j][1] << ")";
            if (j < num_samples - 1) { output << '\t'; }
        }
        output << endl;
    }

    output.close();
    cout << "Total steps: " << res*(num_samples*steps_per_sample + num_exchanges*steps_per_exchange + equilibration_steps) << endl;
}

int main(int argc, char* argv[]) {    
    string filename = argv[1];
    int N = stoi(argv[2]);
    int num_threads = 4;

    auto start = chrono::high_resolution_clock::now();

    func(filename, N, num_threads);

    auto stop = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int microseconds = duration.count();

    cout << "Duration: " << microseconds/1e6 << endl;

}


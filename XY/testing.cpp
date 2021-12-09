#include "SquareXYModel.cpp"
#include <iostream>
#include<chrono>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 8;
    const int L = 1;
    const float J = 1.;
    const float B = 0.;
    const float Bp = 0.;

    int num_temps = 30;
    vector<float> Ts;
    for (int i = 0; i < num_temps; i++) {
        Ts.push_back(float(i)/num_temps*3.);
    }

    cout << float() << endl;


    const int MCStep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);

    int steps_per_exchange = 10*MCStep;
    int num_exchanges = 100;
    int equilibration_steps = 1000*MCStep;

    auto start = chrono::high_resolution_clock::now();

    auto models = parallel_tempering(model, Ts, steps_per_exchange, num_exchanges, equilibration_steps, 4);

    auto stop = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int microseconds = duration.count();

    int nsteps = num_temps*(steps_per_exchange*num_exchanges + equilibration_steps);
    int nsteps_measured = models[0]->nsteps;

    cout << "Actual nsteps: " << nsteps_measured << endl;
    cout << to_string(nsteps) << " steps took " << to_string(float(microseconds)/1000000.) << "s." << endl;
    cout << to_string(float(nsteps)/(float(microseconds)/1000000)) << " steps/second." << endl;

    vector<vector<float>> data(num_temps, vector<float>(2));
    for (int i = 0; i < num_temps; i++) {
        data[i][0] = Ts[i];
        data[i][1] = models[i]->model->energy();
    }

    write_to_file("Data.txt", data);
}

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

    int num_temps = 20;
    vector<float> Ts;
    for (int i = 0; i < num_temps; i++) {
        Ts.push_back(float(i)/num_temps*3.);
    }



    const int MCStep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);

    int steps_per_exchange = 1000*MCStep;
    int num_exchanges = 100;

    auto start = chrono::high_resolution_clock::now();

    vector<SquareXYModel*> models = parallel_tempering(model, Ts, steps_per_exchange, num_exchanges);

    auto stop = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int microseconds = duration.count();

    int nsteps = Ts.size()*steps_per_exchange*num_exchanges;

    cout << to_string(nsteps) << " steps took " << to_string(float(microseconds)/1000000.) << "s." << endl;
    cout << to_string(float(nsteps)/(float(microseconds)/1000000)) << " steps/second." << endl;

    vector<vector<float>> data(num_temps, vector<float>(2));
    for (int i = 0; i < num_temps; i++) {
        data[i][0] = Ts[i];
        data[i][1] = models[i]->energy();
    }

    write_to_file("Data.txt", data);
}

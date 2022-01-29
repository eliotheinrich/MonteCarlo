#include "SquareIsingModel.cpp"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {    
    string filename = argv[1];
    srand((unsigned)time( NULL ));

    const int N = 32;
    const int L = 1;
    const float J = 1.;
    const float B = 0.5;

    const float Tmax = 3.0;
    const float Tmin = 0.1;


    const int MCStep = N*N*L;



    SquareIsingModel *model = new SquareIsingModel(N, L, J, B);
    MonteCarlo<SquareIsingModel> *m = new MonteCarlo<SquareIsingModel>(model);

    int res = 30;
    vector<float> Ms(res);
    vector<float> Es(res);

    int num_samples = 50;
    int steps_per_sample = 50;

    float M;
    float E;
    float T;
    for (int i = 0; i < res; i++) {
        cout << i << endl;
        model->randomize_spins();
        T = float(i)/res*Tmax + float(res - i)/res*Tmin;
        run_MC(m, 5000*MCStep, "trig", T, T, false);

        M = 0.;
        E = 0.;
        for (int j = 0; j < num_samples; j++) {
            M += model->get_magnetization();
            E += model->energy();
            run_MC(m, steps_per_sample*MCStep, "const", T, T, false);
        }
        Ms[i] = M/num_samples;
        Es[i] = model->energy()/num_samples;
    } 

    ofstream output(filename);
    for (int i = 0; i < res; i++ ) {
        output << float(i)/res*Tmax + float(res - i)/res*Tmin << "\t" << Ms[i] << "\t" << Es[i] << endl;
    }
    output.close();

}

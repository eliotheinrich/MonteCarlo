#include "SquareModel.cpp"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {    
    string filename = argv[1];
    srand((unsigned)time( NULL ));

    int N = stoi(argv[2]);
    cout << "N = " << N << endl;
    const int L = 1;
    const float J = 1.;
    const float A = 2.;

    const int MCStep = N*N*L;

    SquareModel *model = new SquareModel(N, L, J, A);
    MonteCarlo<SquareModel> *m = new MonteCarlo<SquareModel>(model);


    const float Tmax = 3.;
    const float Tmin = 0.1;
    int res = 20;

    int num_samples = 5000;
    int steps_per_sample = 5*MCStep;


    vector<float> Ts(res);
    ofstream output(filename);

    for (int i = 0; i < res; i++) {
        Ts[i] = float(i)/res*Tmax + float(res - i)/res*Tmin;
    }

    auto models = parallel_tempering(model, Ts, 100, 2*MCStep, 100*MCStep);

    float T;
    vector<double> sample;

    // Write header
    output << res << "\t" << num_samples << endl;
    for (int i = 0; i < res; i++) {
        T = Ts[i];
        output << Ts[i] << "\t";
        for (int j = 0; j < num_samples; j++) {
            sample = models[i]->model->twist_stiffness();
            models[i]->steps(steps_per_sample, T);
            output << "(" << sample[0] << ", " << sample[1] << ")";
            if (j < num_samples - 1) { output << '\t'; }
        }
        output << endl;
    } 

    output.close();

}

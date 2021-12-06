#include "SquareXYModel.cpp"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {    
    string filename = argv[1];
    srand((unsigned)time( NULL ));

    int N = stoi(argv[2]);
    cout << "N = " << N << endl;
    const int L = 1;
    const float J = 1.;
    const float B = 0.;
    const float Bp = 0.;

    const int MCStep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);
    //TrigonalXYModel *model = new TrigonalXYModel(N, L, J);


    const float Tmax = 3.;
    const float Tmin = 0.1;
    int res = 30;
    vector<float> rhos(res);

    int num_samples = 1000;
    int steps_per_sample = 10;


    float rho;
    float T;
    ofstream output(filename);
    for (int i = 0; i < res; i++) {
        if (i % 5 == 0) {
            cout << i*100/res << "\% finished" << endl;
        }
        model->randomize_spins();
        T = float(i)/res*Tmax + float(res - i)/res*Tmin;
        run_MC(model, 2000*MCStep, "trig", T, T, false);

        rho = 0;
        for (int j = 0; j < num_samples; j++) {
            rho += model->spin_stiffness(T);
            run_MC(model, steps_per_sample*MCStep, "const", T, T, false);
        }

        rhos[i] = rho/num_samples;

        output << T << "\t" << rhos[i] << endl;
    } 

    output.close();


}

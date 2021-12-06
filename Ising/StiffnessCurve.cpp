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


    const float Tmax = 1.;
    const float Tmin = 0.1;
    int res = 30;
    vector<float> rhos(res);

    int num_samples = 50;
    float s1;
    float s2;
    int steps_per_sample = 50;


    float T;
    for (int i = 0; i < res; i++) {
        if (i % 5 == 0) {
            cout << float(i)/res*100 << "\% finished" << endl;
        }
        model->randomize_spins();
        T = float(i)/res*Tmax + float(res - i)/res*Tmin;
        run_MC(model, 5000*MCStep, "trig", T, T, false);

        s1 = 0;
        s2 = 0;
        for (int j = 0; j < num_samples; j++) {
            s1 += model->s1();
            s2 += model->s2();
            run_MC(model, steps_per_sample*MCStep, "const", T, T, false);
        }
        //cout << "s1 = " << s1 << endl;
        //cout << "s2 = " << s2 << endl;
        //cout << "M = " << model->get_magnetization().norm() << endl;
        rhos[i] = J*(s1 - (J/T)*s2)/(num_samples*N*N);
    } 

    ofstream output(filename);
    for (int i = 0; i < res; i++ ) {
        output << float(i)/res*Tmax + float(res - i)/res*Tmin << "\t" << rhos[i] << endl;
    }
    output.close();


}

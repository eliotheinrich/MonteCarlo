#include "SimpleGraphModel.cpp"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {    
    string filename = argv[1];
    srand((unsigned)time( NULL ));

    const int N = 32;
	const float J = 1.;

    const float Tmax = 3.0;
    const float Tmin = 0.1;


    const int MCStep = N;



    SimpleGraphModel *model = new SimpleGraphModel(N, J);
    MonteCarlo<SimpleGraphModel> *m = new MonteCarlo<SimpleGraphModel>(model);

    int res = 30;
    vector<float> Cs(res);
    vector<float> Es(res);

    int num_samples = 50;
    int steps_per_sample = 50;

    float C;
    float E;
    float T;
    for (int i = 0; i < res; i++) {
        cout << i << endl;
        model->init();
        T = float(i)/res*Tmax + float(res - i)/res*Tmin;
        run_MC(m, 5000*MCStep, "trig", T, T, false);

        M = 0.;
        E = 0.;
        for (int j = 0; j < num_samples; j++) {
            C += model->get_connectivity();
            E += model->energy();
            run_MC(m, steps_per_sample*MCStep, "const", T, T, false);
        }
        Cs[i] = C/num_samples;
        Es[i] = model->energy()/num_samples;
    } 

    ofstream output(filename);
    for (int i = 0; i < res; i++ ) {
        output << float(i)/res*Tmax + float(res - i)/res*Tmin << "\t" << Cs[i] << "\t" << Es[i] << endl;
    }
    output.close();

}
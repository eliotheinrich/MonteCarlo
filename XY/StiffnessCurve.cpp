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
    cout << MCStep << endl;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);
    model->random_selection = true;
    MonteCarlo<SquareXYModel> *m = new MonteCarlo<SquareXYModel>(model);
    //TrigonalXYModel *model = new TrigonalXYModel(N, L, J);


    const float Tmax = 3.;
    const float Tmin = 0.1;
    int res = 30;
    vector<float> rhos(res);

    int num_samples = 4000;
    int steps_per_sample = 10;


    float rho;
    vector<float> Ts(res);
    ofstream output(filename);

    for (int i = 0; i < res; i++) {
        Ts[i] = float(i)/res*Tmax + float(res - i)/res*Tmin;
    }

    vector<SquareXYModel*> models = parallel_tempering(model, Ts, 100, 100, 4000*MCStep);
    vector<thread> threads(res);

    function<float(SquareXYModel*)> stiffness;
    function<void(float, int, int)> func;
    float T;
    for (int i = 0; i < res; i++) {
        T = Ts[i];
        stiffness = [T](SquareXYModel *m) { return m->spin_stiffness(T); };
        m = new MonteCarlo<SquareXYModel>(models[i]);
        func = [i, &rhos, m, stiffness](float T, int n1, int n2) { rhos[i] = m->expectation(stiffness, T, n1, n2); };
        threads[i] = thread(func, T, num_samples, steps_per_sample);;
        //rhos[i] = m->expectation(stiffness, T, num_samples, steps_per_sample);
    } 

    for (int i = 0; i < res; i++) {
        threads[i].join(); 
        output << Ts[i] << "\t" << rhos[i] << endl;
    }


    output.close();

}

#include "../TrigonalModel.cpp"
#include <iostream>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 12;
    const int L = 2;
    const float J1 = 1.;
    const float J2 = 0.5;
    const float K1 = 2.0;
    const float K2 = 0.;
    const float K3 = 0.0;
    Vector3f B; B << 0., 0., 0.;
    const float T = 2.0;


    const int MCStep = N*N*L;

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
    //MonteCarlo<TrigonalModel> *m = new MonteCarlo<TrigonalModel>(model);
    //m->steps(10000*MCStep, T);

    for (int i = 0; i < model->V; i++) {
        model->spins[i] << 1., 0., 0.;
    }

    for (int n1 = 0; n1 < N; n1++) {
        model->spins[model->flat_idx(n1, 2, 0, 0)] << 0., 1., 0.;
        model->spins[model->flat_idx(n1, 6, 1, 0)] << 0., 1., 0.;
    }

    int nsteps = 100000;
    int num_frames = 200;
    for (int i = 0; i < nsteps; i++) {
        if (num_frames*i % nsteps == 0) { 
            model->save_spins("data/spins" + to_string(i) + ".txt");
        }
        model->dynamic_step(0.0001);
    }
}

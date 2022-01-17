#include "EasyPlaneHeis.cpp"
#include <iostream>

int main() {
    std::srand((unsigned)std::time( NULL ));

    const int N = 8;
    const int L = 1;
    const float J = 1.;
    const float K = 0.3;
    const float T = 0.1;


    const int MCStep = N*N*L;

    EasyPlaneHeis *model = new EasyPlaneHeis(N, L, J, 0.);

    MonteCarlo<EasyPlaneHeis> *m = new MonteCarlo<EasyPlaneHeis>(model);

    int num_samples = 1000;
    std::vector<std::vector<float>> log = std::vector<std::vector<float>>(num_samples, std::vector<float>(2));
    for (int i = 0; i < num_samples; i++) {
        m->steps(MCStep, T);
        log[i][0] = model->get_magnetization().norm();
        log[i][1] = model->energy();
    }

    std::ofstream output("Log.txt");
    for (int i = 0; i < log.size(); i++) {
        output << log[i][0] << "\t" << log[i][1] << "\n";
    }

    model->save_spins("Spins.txt");

    std::cout << "Energy = " << model->energy() << std::endl;


}

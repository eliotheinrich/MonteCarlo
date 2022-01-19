#include "TrigonalModel.cpp"
#include <iostream>

int main() {
    std::srand((unsigned)std::time( NULL ));

    const int N = 8;
    const int L = 1;
    const float J1 = 1.;
    const float J2 = 0.3;
    const float K1 = 2.;
    const float K2 = 0.;
    const float K3 = 0.15;
    Eigen::Vector3f B; B << 0., 0., 0.;
    const float T = .1;


    int MCStep;
    bool cluster = true;
    if (cluster) {
        MCStep = 1;
    } else {
        MCStep = N*N*L;
    }

    TrigonalModel *model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
    MonteCarlo<TrigonalModel> *m = new MonteCarlo<TrigonalModel>(model);

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

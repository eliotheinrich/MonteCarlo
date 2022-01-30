#include "TrigonalModel.cpp"
#include <iostream>

class TrigonalModelM : public TrigonalModel {
    public:
        TrigonalModelM(int N, int L, float J1, float J2, float K1, float K2, float K3,
                                    Eigen::Vector3d B) : TrigonalModel(N, L, J1, J2, K1, K2, K3, B) {}

        TrigonalModelM* clone() {
            TrigonalModelM* new_model = new TrigonalModelM(N, L, J1, J2, K1, K2, K3, B);
            return new_model;
        }

        std::vector<double> tracking_func(int i) {
            Eigen::Vector3d M = spins[i]/V;
            return std::vector<double>{
                onsite_energy(i) + 2*bond_energy(i),
                M[0], M[1], M[2]
            };
        }

        std::vector<double> init_func() {
            Eigen::Vector3d M = get_magnetization();
            return std::vector<double>{
                energy(),
                M[0], M[1], M[2]
            };
        }
};

int main() {
    std::srand((unsigned)std::time( NULL ));

    const int N = 16;
    const int L = 1;
    const float J1 = .64;
    const float J2 = 0.;
    const float K1 = 2.56;
    const float K2 = 0.;
    const float K3 = 0.15;
    Eigen::Vector3d B; B << 0., 0., 0.;
    const float T = .4;


#ifdef CLUSTER_UPDATE
    const int MCStep = 1;
#else
    const int MCStep = N*N*L;
#endif

    TrigonalModelM *model = new TrigonalModelM(N, L, J1, J2, K1, K2, K3, B);
    MonteCarlo<TrigonalModelM> *m = new MonteCarlo<TrigonalModelM>(model);

    model->start_tracking();
    int num_samples = 4000*MCStep;
    std::vector<std::vector<float>> log = std::vector<std::vector<float>>(num_samples, std::vector<float>(2));
    for (int i = 0; i < num_samples; i++) {
        m->steps(1, T);
        log[i][0] = sqrt(pow(model->q[1],2) + pow(model->q[2],2) + pow(model->q[3],2));
        log[i][1] = model->q[0];
    }

    std::ofstream output("Log.txt");
    for (int i = 0; i < log.size(); i++) {
        output << log[i][0] << "\t" << log[i][1] << "\n";
    }

    model->save_spins("Spins.txt");

    std::cout << "Energy = " << m->model->energy() << std::endl;
    std::cout << "Magnetization = " << m->model->get_magnetization().norm() << std::endl;
}

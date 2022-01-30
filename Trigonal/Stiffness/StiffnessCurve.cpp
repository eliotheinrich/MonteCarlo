#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

class TrigonalModelS : public TrigonalModel {
    public:
        TrigonalModelS(int N, int L, float J1, float J2, float K1, float K2, float K3,
                                    Eigen::Vector3d B) : TrigonalModel(N, L, J1, J2, K1, K2, K3, B) {}

        TrigonalModelS* clone() {
            TrigonalModelS* new_model = new TrigonalModelS(N, L, J1, J2, K1, K2, K3, B);
            return new_model;
        }


        std::vector<double> tracking_func(int i) {
            auto twist = twist_derivatives(i);
            for (int j = 0; j < twist.size(); j++) {
                twist[j] *= 2.;
            }
            return twist;
        }

        std::vector<double> init_func() {
            return twist_derivatives();
        }
};

std::vector<double> twisting_sampler(TrigonalModelS *model) {
    return SpinModel::twist_terms(model->q);
}

int main(int argc, char* argv[]) {
    srand((unsigned)time( NULL ));

    std::string filename = argv[1];
    int num_threads = std::stoi(argv[2]);
    int N = std::stoi(argv[3]);
    int L = std::stoi(argv[4]);

    std::cout << "N = " << N << std::endl;
    std::cout << "L = " << L << std::endl;

    // Base case
    float J1 = std::stof(argv[5]);
    float J2 = std::stof(argv[6]);
    float K1 = std::stof(argv[7]);
    float K2 = std::stof(argv[8]);
    float K3 = std::stof(argv[9]);
    Eigen::Vector3d B; B << std::stof(argv[10]), std::stof(argv[11]), std::stof(argv[12]);

#ifdef CLUSTER_UPDATE
    const int MCStep = 1;
#else
    const int MCStep = N*N*L;
#endif


    TrigonalModelS *model;
    /*
    MonteCarlo<TrigonalModel> *m1;
    for (int i = 0; i < 5; i++) {
        model = new TrigonalModel(N, L, J1, J2, K1, K2, K3, B);
        m1 = new MonteCarlo<TrigonalModel>(model);
        std::cout << "Pre: " << m1->model->energy() << std::endl;
        m1->steps(MCStep*10000, 0.5);
        std::cout << "Post: " << m1->model->energy() << std::endl;
    }*/

    model = new TrigonalModelS(N, L, J1, J2, K1, K2, K3, B);
    /*
    MonteCarlo<TrigonalModel> *m2 = new MonteCarlo<TrigonalModel>(model->clone());
    std::cout << "Post: " << m1->model->energy() << std::endl;
    std::cout << "Pre: " << m2->model->energy() << std::endl;
    m2->steps(MCStep*10000, 0.5);
    std::cout << "Post: " << m2->model->energy() << std::endl;
    */


    const float Tmax = 2.;
    const float Tmin = 0.1;
    unsigned int resolution = std::stoi(argv[13]);

    std::vector<double> T(resolution);

    for (int i = 0; i < resolution; i++) {
        T[i] = double(i)/resolution*Tmax + double(resolution - i)/resolution*Tmin;
    }

    unsigned long long steps_per_run = std::stoi(argv[14])*MCStep;

    unsigned int num_samples = std::stoi(argv[15]);
    unsigned long long steps_per_sample = std::stoi(argv[16])*MCStep;
    steps_per_sample = N*N*L*10;

    int num_runs = std::stoi(argv[17]);

    unsigned long long nsteps = resolution*num_runs*(steps_per_run + num_samples*steps_per_sample);
    std::cout << "nsteps: " << nsteps << std::endl;
    std::cout << "Expected completion time: " << (long long) 2*nsteps/3300000./num_threads/60. << " minutes. " << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto data = sample_r(twisting_sampler, model, T, steps_per_run, num_samples, steps_per_sample, num_threads);
    std::string header = std::to_string(N) + "\t" + std::to_string(L);
    write_data(&data, &T, filename, header);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count()/1000000.;

    std::cout << "Completion time: " << seconds/60. << " minutes." << std::endl;


}


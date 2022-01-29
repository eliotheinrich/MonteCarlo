#include "../SquareXYModel.cpp"
#include "../../Utility.cpp"
#include "../../Routines.cpp"
#include <iostream>
#include <math.h>

int main(int argc, char* argv[]) {
    std::srand((unsigned)std::time( NULL ));

    int N = 32;
    const int L = 1;
    const float J = 1.;
    const float T = 1.5;

#ifdef CLUSTER_UPDATE
    const int MCStep = 1;
#else
    const int MCStep = N*N*L;
#endif
    std::cout << "N = " << N << std::endl;

    SquareXYModel *model = new SquareXYModel(N, L, J, 0., 0.);
    MonteCarlo<SquareXYModel> *m =  new MonteCarlo<SquareXYModel>(model);

    m->steps(2000*MCStep, T);

    std::vector<float> energy_density(model->V);
    for (int i = 0; i < model->V; i++) {
        energy_density[i] = model->onsite_energy(i) + model->bond_energy(i);
    }

    std::ofstream output("density.txt");
    output << N << "\n";
    for (int i = 0; i < model->V; i++) {
        output << energy_density[i] << "\t";
    }
    output << std::endl;

    std::cout << model->energy() << std::endl;


}


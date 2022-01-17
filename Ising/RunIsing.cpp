#include "SquareIsingModel.cpp"
#include <iostream>

int main() {
    srand((unsigned)time( NULL ));

    const int N = 20;
    const int L = 1;
    const float J = 1.;
    const float B = 0.;
    const float T = 0.1;


    const int MCStep = N*N*L;

    SquareIsingModel *model = new SquareIsingModel(N, L, J, B);
    MonteCarlo<SquareIsingModel> *m = new MonteCarlo<SquareIsingModel>(model);

    std::cout << model->get_magnetization() << std::endl;
    m->steps(50000*MCStep, T);
    std::cout << model->get_magnetization() << std::endl;

    model->save_spins("Spins.txt");
    

}

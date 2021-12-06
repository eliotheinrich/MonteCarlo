#include "SquareXYModel.cpp"
#include <iostream>

using namespace std;

int main() {
    srand((unsigned)time( NULL ));

    const int N = 6;
    const int L = 1;
    const float J = 1.;
    const float B = 0.;
    const float Bp = 0.;
    const float T = 0.01;


    const int MCStep = N*N*L;

    SquareXYModel *model = new SquareXYModel(N, L, J, B, Bp);
    for (int n1 = 0; n1 < N; n1++) {
        for (int n2 = 0; n2 < N; n2++) {
            for (int n3 = 0; n3 < L; n3++) {
                model->spins[n1][n2][n3][0] << 1.,0.;
            }
        }
    }

    float E1 = model->energy();

    model->spins[N/2][N/2][0][0] << 0.,1.;
    model->spins[N/2+1][N/2][0][0] << -1.,0.;
    model->spins[N/2+1][N/2+1][0][0] << 0.,-1.;

    run_MC(model, 10*N*N*L, "const", T, T);

    float E2 = model->energy();
    
    cout << "Chemical potential: " << E2 - E1 << endl;

    model->save_spins("Spins.txt");
    

}

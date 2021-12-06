#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../XYModel.cpp"
#include "../Utility.cpp"


using namespace std;
using namespace Eigen;

class TrigonalXYModel : public XYModel {
    public:
        int N;
        int L;
        float J;

        TrigonalXYModel(int N, int L, float J) : XYModel(1, N, N, L) {
            this->N = N;
            this->L = L;
            this->J = J;

            function<float(Vector2f, Vector2f)> dotfunc = [J](Vector2f S1, Vector2f S2) {
                return -J*S1.dot(S2);
            };

            this->add_bond(Bond{1,0,0,0,dotfunc});
            this->add_bond(Bond{-1,0,0,0,dotfunc});
            this->add_bond(Bond{0,1,0,0,dotfunc});
            this->add_bond(Bond{0,-1,0,0,dotfunc});
            this->add_bond(Bond{1,-1,0,0,dotfunc});
            this->add_bond(Bond{-1,1,0,0,dotfunc});
        }

        const float onsite_energy(int n1, int n2, int n3, int s) {
            return 0.;
        }

        const float bond_energy(int n1, int n2, int n3, int s) {
            float E = 0;

            // NN interactions

            // PBC boundary conditions
            E -= this->spins[n1][n2][n3][s].dot(this->spins[n1][(n2+1)%N][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[n1][mod(n2-1, N)][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[(n1+1)%N][n2][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[mod(n1-1, N)][n2][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[(n1+1)%N][mod(n2-1, N)][n3][s]);
            E -= this->spins[n1][n2][n3][s].dot(this->spins[mod(n1-1, N)][(n2+1)%N][n3][s]);

            return 0.5*J*E;
        }
};


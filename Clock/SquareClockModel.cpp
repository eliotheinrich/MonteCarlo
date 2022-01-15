#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include <complex>
#include "../ClockModel.cpp"
#include "../Utility.cpp"


using namespace std;
using namespace Eigen;

template <int q>
class SquareClockModel : public ClockModel<q> {
    public:
        int N;
        int L;
        float J;

        vector<vector<float>> bond_table;

        SquareClockModel(int N, int L, float J) : ClockModel<q>(N, N, L) {
            this->N = N;
            this->L = L;
            this->J = J;


            this->bond_table = vector<vector<float>>(q, vector<float>(q));

            for (int i = 0; i < q; i++) {
                for (int j = 0; j < q; j++) {
                    bond_table[i][j] = -J*cos(2*PI/q*(i - j));
                }
            }

            function<float(int, int)> bondfunc = [J, this](int p1, int p2) {
                return this->bond_table[p1][p2];
            };



            Vector3f v1; v1 << 1.,0.,0.;
            Vector3f v2; v2 << 0.,1.,0.;
            this->add_bond(1,0,0,   v1, bondfunc);
            this->add_bond(-1,0,0, -v1, bondfunc);
            this->add_bond(0,1,0,   v2, bondfunc);
            this->add_bond(0,-1,0, -v2, bondfunc);
        }

        SquareClockModel* clone() {
            SquareClockModel<q>* new_model = new SquareClockModel<q>(N, L, J);
            for (int i = 0; i < this->V; i++) {
                new_model->spins[i]  = this->spins[i];
            }
            return new_model;
        }

        const float onsite_energy(int i) {
            return 0.;
        }
};


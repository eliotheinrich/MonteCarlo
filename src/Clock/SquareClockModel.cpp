#ifndef SQUARECLOCK_
#define SQUARECLOCK_

#include "SquareClockModel.h"
#include <iostream>
#include <functional>
#include <vector>
#include <string>

template <int q>
SquareClockModel<q>::SquareClockModel(Params &params) {
    this->N = params.geti("system_size");
    this->L = params.geti("layers", DEFAULT_LAYERS);
    ClockModel<q>::init_params(N, N, L);

    this->J = params.getf("J");

    for (int i = 0; i < q; i++) {
        for (int j = 0; j < q; j++) {
            bond_table[i][j] = -J*cos(2*PI/q*(i - j));
        }
    }

    std::function<float(int, int)> bondfunc = [this](int p1, int p2) {
        return this->bond_table[p1][p2];
    };



    Eigen::Vector3d v1; v1 << 1.,0.,0.;
    Eigen::Vector3d v2; v2 << 0.,1.,0.;
    this->add_bond(1,0,0,   v1, bondfunc);
    this->add_bond(-1,0,0, -v1, bondfunc);
    this->add_bond(0,1,0,   v2, bondfunc);
    this->add_bond(0,-1,0, -v2, bondfunc);
}

template <int q>
double SquareClockModel<q>::onsite_energy(int i) const {
    return 0.;
}


#endif
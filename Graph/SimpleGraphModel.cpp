#ifndef SIMPLEGRAPH_
#define SIMPLEGRAPH_

#include "SimpleGraphModel.h"

SimpleGraphModel::SimpleGraphModel(int N, float J) : GraphModel(N) {
    this->N = N;
    this->J = J;
}

SimpleGraphModel* SimpleGraphModel::clone() {
    SquareIsingModel *new_model = new SimpleGraphModel(N, J);
}

const float SimpleGraphModel::onsite_energy(int i) {
    float E = 0;

    return E;
}

const float SimpleGraphModel::bond_energy(int i) {
    float E = 0;
    E += deg(i);
    return J*E;
}

#endif
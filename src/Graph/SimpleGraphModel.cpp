#include "SimpleGraphModel.h"

SimpleGraphModel::SimpleGraphModel(Params &params) {
    this->N = get<int>(params, "system_size");
    GraphModel::init_params(N);
    this->J = get<double>(params, "J");
}

double SimpleGraphModel::onsite_energy(int i) const {
    double E = 0;

    return E;
}

double SimpleGraphModel::bond_energy(int i) const {
    double E = 0;
    E += deg(i);
    return J*E;
}
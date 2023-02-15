#include "SquareIsingModel.h"

SquareIsingModel::SquareIsingModel(int N, int L, float J, float B) : IsingModel(N, N, L) {
    this->N = N;
    this->L = L;
    this->J = J;
    this->B = B;
}

SquareIsingModel* SquareIsingModel::clone() {
    SquareIsingModel *new_model = new SquareIsingModel(N, L, J, B);
    for (int i = 0; i < V; i++) {
        new_model->spins[i] = this->spins[i];
    }
    return new_model;
}

double SquareIsingModel::onsite_energy(int i) const {
    float E = 0;

    // Onsite interactions
    E -= B*this->spins[i];

    return E;
}

double SquareIsingModel::bond_energy(int i) const {
    float E = 0;

    Eigen::Vector3i idxs = tensor_idx(i);
    int n1 = idxs[0]; int n2 = idxs[1]; int n3 = idxs[2]; 

    // NN interactions
    E -= J*spins[i]*spins[flat_idx(mod(n1+1, N), n2, n3)];
    E -= J*spins[i]*spins[flat_idx(mod(n1-1, N), n2, n3)];

    E -= J*spins[i]*spins[flat_idx(n1, mod(n2+1, N), n3)];
    E -= J*spins[i]*spins[flat_idx(n1, mod(n2-1, N), n3)];

    return 0.5*E;
}

std::map<std::string, int> SquareIsingModel::get_int_params() const {
    std::map<std::string, int> data;
    data.emplace("system_size", this->N);
    data.emplace("layers", this->L);
    return data;
}

std::map<std::string, double> SquareIsingModel::get_double_params() const {
    std::map<std::string, double> data;
    data.emplace("T", this->T);
    data.emplace("J", this->J);
    data.emplace("B", this->B);
    return data;
}
#ifndef ISINGMC_
#define ISINGMC_

#include "IsingModel.h"

IsingModel::IsingModel(int N1, int N2, int N3) {
    this->N1 = N1;
    this->N2 = N2;
    this->N3 = N3;
    this->V = N1*N2*N3;

    this->spins = std::vector<float>(V);

    this->randomize_spins();

    this->acceptance = 0.5;
}

inline int IsingModel::flat_idx(int n1, int n2, int n3) const {
    return n1 + N1*(n2 + N2*n3);
}

inline Eigen::Vector3i IsingModel::tensor_idx(int i) const {
    int n1 = i % N1;
    i = i / N1;
    int n2 = i % N2;
    i = i / N2;
    int n3 = i % N3;
    Eigen::Vector3i v; v << n1, n2, n3;
    return v;
}

void IsingModel::randomize_spins() {
    for (int i = 0; i < V; i++) {
        // For each site, initialize spin randomly
        if (rand() % 2) {
            this->spins[i] = 1.;
        } else {
            this->spins[i] = -1.;
        }
    }
}

inline double IsingModel::get_magnetization() const {
    double M = 0;
    for (int i = 0; i < V; i++) {
        M += this->spins[i];
    }
    
    return M/(N1*N2*N3);
}

void IsingModel::generate_mutation() {
    // Randomly select a site to mutate
    mut.i = std::rand() % V;
}

void IsingModel::accept_mutation() {
    return;
}

void IsingModel::reject_mutation() {
    this->spins[mut.i] = -this->spins[mut.i];
}

double IsingModel::energy() const {
    float E = 0;

    for (int i = 0; i < V; i++) {
        E += onsite_energy(i);
        E += bond_energy(i);
    }
    return E;
}

double IsingModel::energy_change() {
    float E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
    this->spins[mut.i] = -this->spins[mut.i];
    float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

    return E2 - E1;
}

// Saves current spin configuration
void IsingModel::save_spins(std::string filename) {
    std::ofstream output_file;
    output_file.open(filename);
    output_file << N1 << "\t" << N2 << "\t" << N3 << "\n";
    for (int i = 0; i < V; i++) {
        output_file << spins[i];
        if (i < V-1) { output_file << "\t"; }
    }
    output_file.close();
}        

#endif

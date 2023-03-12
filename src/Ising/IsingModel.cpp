#include "IsingModel.h"
#include <fstream>
#include <math.h>

void IsingModel::init_params(int N1, int N2, int N3) {
    this->N1 = N1;
    this->N2 = N2;
    this->N3 = N3;
    this->V = N1*N2*N3;

    this->acceptance = 0.5;
}

void IsingModel::init() {
    this->spins = std::vector<float>(V);
    this->randomize_spins();
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

double IsingModel::get_magnetization() const {
    double M = 0;
    for (int i = 0; i < V; i++) {
        M += this->spins[i];
    }
    
    return M/(N1*N2*N3);
}

void IsingModel::generate_mutation() {
    // Randomly select a site to mutate
    mut.i = rand() % V;
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

std::map<std::string, Sample> IsingModel::take_samples() const {
    std::map<std::string, Sample> samples;
    samples.emplace("energy", energy());
    samples.emplace("magnetization", get_magnetization());
    return samples;
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
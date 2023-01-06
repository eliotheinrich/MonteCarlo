#ifndef CLOCKMC_
#define CLOCKMC_

#include "ClockModel.h"

template <int q>
ClockModel<q>::ClockModel(int N1, int N2, int N3) {
    this->N1 = N1;
    this->N2 = N2;
    this->N3 = N3;
    this->V = N1*N2*N3;

    this->spins = std::vector<int>(V);
    this->neighbors = std::vector<std::vector<int>>(V, std::vector<int>(0));

    this->r.seed(rand());
    this->randomize_spins();

    this->mut.i = 0;
    this->mut_mode = 0;
}


template <int q>
inline const int ClockModel<q>::flat_idx(int n1, int n2, int n3) {
    return n1 + N1*(n2 + N2*n3);
}

template <int q>
inline const Eigen::Vector3i ClockModel<q>::tensor_idx(int i) {
    int n1 = i % N1;
    i = i / N1;
    int n2 = i % N2;
    i = i / N2;
    int n3 = i % N3;
    Eigen::Vector3i v; v << n1, n2, n3;
    return v;
}

template <int q>
void ClockModel<q>::randomize_spins() {
    for (int i = 0; i < V; i++) {
        // For each site, initialize spin randomly
        spins[i] = r() % q;
    }
}

template <int q>
void ClockModel<q>::add_bond(int d1, int d2, int d3, Eigen::Vector3d v, std::function<float(int, int)> bondfunc) {
    ClockBond b{d1, d2, d3, v, bondfunc};
    this->bonds.push_back(b);
    int i; int j;
    for (int n1 = 0; n1 < N1; n1++) {
        for (int n2 = 0; n2 < N2; n2++) {
            for (int n3 = 0; n3 < N3; n3++) {
                i = flat_idx(n1, n2, n3);
                j = flat_idx(mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3));
                neighbors[i].push_back(j);
            }
        }
    }
}

template <int q>
inline float ClockModel<q>::get_magnetization() {
    float M = 0;
    float x = 0; float y = 0;
    for (int i = 0; i < V; i++) {
        x += cos(2*PI*spins[i]/q);
        y += sin(2*PI*spins[i]/q);
    }
    
    return sqrt(x*x + y*y)/(N1*N2*N3);
}

template <int q>
void ClockModel<q>::metropolis_mutation() {
    mut.dq = r() % 3 - 1;
}

template <int q>
void ClockModel<q>::cluster_update() {
    s.clear();

    std::stack<int> c;
    int m = r() % V;
    c.push(m);

    int p = r() % q;

    int j; float dE;
    int s_new; 
    while (!c.empty()) {
        m = c.top();
        c.pop();

        if (!s.count(m)) {
            s.insert(m);

            s_new = mod(2*spins[m] - p, q);

            for (int n = 0; n < neighbors[m].size(); n++) {
                j = neighbors[m][n];

                if (!s.count(j)) {
                    dE = bonds[n].bondfunc(spins[j], s_new) - bonds[n].bondfunc(spins[j], spins[m]);
                    if ((float) r()/RAND_MAX < 1. - std::exp(-dE/T)) {
                        c.push(j);
                    }
                }
            }
            spins[m] = s_new;
        }
    }
}

template <int q>
void ClockModel<q>::generate_mutation() {
    mut.i++;
    if (mut.i == V) {
        mut.i = 0;
        mut_mode++;
    }

    if (mut_mode < 3) {
        metropolis_mutation();
    } else {
        cluster_update();
        mut_mode = 0;
    }
}

template <int q>
void ClockModel<q>::accept_mutation() {
    return;
}

template <int q>
void ClockModel<q>::reject_mutation() {
    this->spins[mut.i] = mod(this->spins[mut.i] - mut.dq, q);
}

template <int q>
const float ClockModel<q>::bond_energy(int i) {
    float E = 0.;
    int j;
    for (int n = 0; n < bonds.size(); n++) {
        j = neighbors[i][n];
        E += 0.5*bonds[n].bondfunc(spins[i], spins[j]);
    }

    return E;
}

template <int q>
const float ClockModel<q>::energy() {
    float E = 0;

    for (int i = 0; i < V; i++) {
        E += onsite_energy(i);
        E += bond_energy(i);
    }
    return E;
}

template <int q>
const float ClockModel<q>::energy_change() {
    if (mut_mode == 0) { 
        return -1.;
    }
    else {
        float E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
        this->spins[mut.i] = mod(this->spins[mut.i] + mut.dq, q);
        float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
        return E2 - E1;
    }
}

// Saves current spin configuration
template <int q>
void ClockModel<q>::save_spins(std::string filename) {
    std::ofstream output_file(filename);

    output_file << N1 << "\t" << N2 << "\t" << N3 << "\n";
    for (int i = 0; i < V; i++) {
        output_file << cos(2*PI*spins[i]/q) << "\t" << sin(2*PI*spins[i]/q);
        if (i < V - 1) { output_file << "\t"; }
    }
    output_file.close();
}        

#endif

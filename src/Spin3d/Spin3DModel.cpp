#include "Spin3DModel.h"
#include <stack>
#include <iostream>
#include <fstream>

GaussianDist::GaussianDist(float mean, float std) {
	this->rd.seed(rand());
	this->gen = std::default_random_engine(rd());
	this->dist = std::normal_distribution<>(mean, std);
}

float GaussianDist::sample() {
	return dist(gen);
}

void Spin3DModel::init_params(int sl, int N1, int N2=-1, int N3=-1) {
    this->sl = sl;
    this->N1 = N1;
    if (N2 == -1) { this->N2 = N1; } else { this->N2 = N2; }
    if (N3 == -1) { this->N3 = N1; } else { this->N3 = N3; }
    this->V = N1*N2*N3*sl;

    this->spins = std::vector<Eigen::Vector3d>(V);
#ifdef CLUSTER_UPDATE
    this->neighbors = std::vector<std::vector<int>>(V+1, std::vector<int>(0));
#else
    this->neighbors = std::vector<std::vector<int>>(V, std::vector<int>(0));
#endif
}

void Spin3DModel::init() {
#ifdef CLUSTER_UPDATE
    // Connect every site to the ghost 
    for (int i = 0; i < V; i++) {
        neighbors[V].push_back(i);
        neighbors[i].push_back(V);
    }

    this->s0 = Eigen::Matrix3d::Identity();
#endif

    this->randomize_spins();

    this->acceptance = 0.5;
    this->sigma = 0.25;

    this->dist = GaussianDist(0., 1.0);

    this->mut.i = 0;
    this->tracking = false;
}

std::vector<double> Spin3DModel::tracking_func(int i) {
    return std::vector<double>(0);
}

std::vector<double> Spin3DModel::init_func() {
    std::vector<double> v = tracking_func(0);
    int dtype_size = v.size();

    v = std::vector<double>(dtype_size, 0);
    std::vector<double> vt(dtype_size);
    for (int i = 0; i < V; i++) {
        vt = tracking_func(i);
        for (int j = 0; j < dtype_size; j++) {
            v[j] += vt[j];
        }
    }

    return v;
}

void Spin3DModel::start_tracking() {
    this->tracking = true;
    this->q = init_func();
}

inline void Spin3DModel::set_spin(int i, Eigen::Vector3d S) {
    if (tracking) {
        std::vector<double> q1 = tracking_func(i);

        spins[i] = S;

        std::vector<double> q2 = tracking_func(i);

        // Update tracked quantities
        for (int j = 0; j < q.size(); j++) {
            q[j] += q2[j] - q1[j];
        }
    } else { 
        spins[i] = S;
    }
}

void Spin3DModel::randomize_spins() {
    for (int i = 0; i < V; i++) {
        spins[i] = Eigen::Vector3d::Random(3).normalized();
    }
}

void Spin3DModel::add_bond(int d1, int d2, int d3, int ds, Eigen::Vector3d v, std::function<double(Eigen::Vector3d, Eigen::Vector3d)> bondfunc) {
    HeisBond b{d1, d2, d3, ds, v, bondfunc};
    this->bonds.push_back(b);
    int i; int j;
    for (int n1 = 0; n1 < N1; n1++) {
        for (int n2 = 0; n2 < N2; n2++) {
            for (int n3 = 0; n3 < N3; n3++) {
                for (int s = 0; s < sl; s++) {
                    i = flat_idx(n1, n2, n3, s);
                    j = flat_idx(mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3), mod(s + b.ds, sl));
                    neighbors[i].push_back(j);
                    std::rotate(neighbors[i].rbegin(), neighbors[i].rbegin()+1, neighbors[i].rend());
                }
            }
        }
    }

    double f = v[0]*alpha;
    Eigen::Matrix3d R;
    R << cos(f), -sin(f), 0,
            sin(f), cos(f), 0.,
            0., 0., 1.;
    R1s.push_back(R);
    R2s.push_back(R*R);
    R3s.push_back(R*R*R);
}

std::vector<double> Spin3DModel::twist_terms(std::vector<double> dE) {
    return std::vector<double>{dE[3], dE[2], dE[0], pow(dE[1],2), dE[1], dE[0]*dE[2], std::pow(dE[0],2)*dE[1], 
                                dE[0]*dE[1], std::pow(dE[0],2), std::pow(dE[0],4), std::pow(dE[0],3)};
}

std::vector<double> Spin3DModel::twist_derivatives(int i) const {
    double E0 = 0.;
    double E1 = 0.;
    double E2 = 0.;
    double E3 = 0.;
    double Em1 = 0.;
    double Em2 = 0.;
    double Em3 = 0.;

    int j;
    Eigen::Vector3d S1 = spins[i];
    Eigen::Vector3d S2;
    for (int n = 0; n < bonds.size(); n++) {
        j = neighbors[i][n];
        S2 = spins[j];

        E0 += bonds[n].bondfunc(S1, S2);

        E1 += bonds[n].bondfunc(S1, R1s[n]*S2);
        Em1 += bonds[n].bondfunc(S1, R1s[n].transpose()*S2);

        E2 += bonds[n].bondfunc(S1, R2s[n]*S2);
        Em2 += bonds[n].bondfunc(S1, R2s[n].transpose()*S2);

        E3 += bonds[n].bondfunc(S1, R3s[n]*S2);
        Em3 += bonds[n].bondfunc(S1, R3s[n].transpose()*S2);
    }

    double d1E = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha/2.;
    double d2E = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/pow(alpha, 2)/2.;
    double d3E = (1./8.*Em3 - 1.*Em2 + 13./8.*Em1 - 13./8.*E1 + 1.*E2 - 1./8.*E3)/pow(alpha, 3)/2.;
    double d4E = (-1./6.*Em3 + 2.*Em2 - 13./2.*Em1 + 28./3.*E0 - 13./2.*E1 + 2.*E2 - 1./6.*E3)/pow(alpha, 4)/2.;


    return std::vector<double>{d1E, d2E, d3E, d4E};
}

std::vector<double> Spin3DModel::twist_derivatives() const {
    std::vector<double> twist = std::vector<double>(4, 0);
    std::vector<double> twist_i = std::vector<double>(4, 0);

    for (int i = 0; i < V; i++) {
        twist_i = twist_derivatives(i);
        for (int j = 0; j < 4; j++) {
            twist[j] += twist_i[j];
        }
    }
    
    return twist;
}

Eigen::Vector3d Spin3DModel::get_magnetization() const {
    Eigen::Vector3d M = Eigen::Vector3d::Constant(0);
    for (int i = 0; i < V; i++) {
        M += spins[i];
    }
    
#ifdef CLUSTER_UPDATE
    return s0.transpose()*M/V;
#else
    return M/V;
#endif
}

std::vector<double> Spin3DModel::correlation_function(int i, int a = 2, int b = 2) const {
    std::vector<double> Cij = std::vector<double>(V); 

    int j;
    Eigen::Vector4i idxs = tensor_idx(i);
    int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];
    for (int n1 = 0; n1 < N1; n1++) {
        for (int n2 = 0; n2 < N2; n2++) {
            for (int n3 = 0; n3 < N3; n3++) {
                for (int s = 0; s < sl; s++) {
                    j = flat_idx(n1, n2, n3, s);
                    Cij[j] = spins[flat_idx((m1 + n1)%N1, 
                                            (m2 + n2)%N2, 
                                            (m3 + n3)%N3, 
                                            (s + k)%sl)][a]*spins[i][b];
                }
            }
        }
    }
    return Cij;
}

std::vector<double> Spin3DModel::full_correlation_function(int i) const {
    std::vector<double> Cij = std::vector<double>(V); 

    int j;
    Eigen::Vector4i idxs = tensor_idx(i);
    int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];
    for (int n1 = 0; n1 < N1; n1++) {
        for (int n2 = 0; n2 < N2; n2++) {
            for (int n3 = 0; n3 < N3; n3++) {
                for (int s = 0; s < sl; s++) {
                    j = flat_idx(n1, n2, n3, s);
                    Cij[j] = spins[flat_idx((m1 + n1)%N1, 
                                            (m2 + n2)%N2, 
                                            (m3 + n3)%N3, 
                                            (s + k)%sl)].dot(spins[i]);
                }
            }
        }
    }
    return Cij;
}

double Spin3DModel::skyrmion_density(int i) const {
    int j;

    Eigen::Vector3d dSdX; dSdX << 0., 0., 0.;
    Eigen::Vector3d dSdY; dSdY << 0., 0., 0.;
    for (int n = 0; n < bonds.size(); n++) {
        j = neighbors[i][n];
        if (bonds[n].v[0] != 0.) {
            dSdX += bonds[n].v[0]*(spins[j] - spins[i]);
        }
        if (bonds[n].v[1] != 0.) {
            dSdX += bonds[n].v[1]*(spins[j] - spins[i]);
        }

    }
    dSdX = dSdX/bonds.size();
    dSdY = dSdY/bonds.size();

    return spins[i].dot(dSdX.cross(dSdY));
}

std::vector<double> Spin3DModel::skyrmion_correlation_function(int i) const {
    std::vector<double> Cij = std::vector<double>(V); 

    int j;

    Eigen::Vector4i idxs = tensor_idx(i);
    int m1 = idxs[0]; int m2 = idxs[1]; int m3 = idxs[2]; int k = idxs[3];

    double Si = skyrmion_density(i);
    for (int n1 = 0; n1 < N1; n1++) {
        for (int n2 = 0; n2 < N2; n2++) {
            for (int n3 = 0; n3 < N3; n3++) {
                for (int s = 0; s < sl; s++) {
                    j = flat_idx(n1, n2, n3, s);
                    Cij[j] = skyrmion_density(flat_idx((m1 + n1)%N1, 
                                                        (m2 + n2)%N2, 
                                                        (m3 + n3)%N3, 
                                                        (s + k)%sl))*Si;
                }
            }
        }
    }
    return Cij;
}

#ifdef CLUSTER_UPDATE
void Spin3DModel::cluster_update() {
    s.clear();

    std::stack<int> c;
    int m = rand() % V;
    c.push(m);

    Eigen::Vector3d ax; ax << dist.sample(), dist.sample(), dist.sample();
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() - 2*ax*ax.transpose()/std::pow(ax.norm(), 2);

    int j; double dE;
    Eigen::Matrix3d s0_new;
    Eigen::Vector3d s_new;
    bool is_ghost; bool neighbor_is_ghost;
    while (!c.empty()) {
        m = c.top();
        c.pop();


        if (!s.count(m)) {
            s.insert(m);

            is_ghost = (m == V);
            if (is_ghost) { // Site is ghost
                s0_new = R*s0;
            } else {
                s_new = R*spins[m];
            }

            for (int n = 0; n < neighbors[m].size(); n++) {
                j = neighbors[m][n];
                neighbor_is_ghost = (j == V);

                if (!s.count(j)) {
                    if (neighbor_is_ghost) {
                        dE = onsite_func(s0.transpose()*s_new) - onsite_func(s0.transpose()*spins[m]);
                    } else if (is_ghost) {
                        dE = onsite_func(s0_new.transpose()*spins[j]) - onsite_func(s0.transpose()*spins[j]);
                    } else { // Normal bond
                        dE = bonds[n].bondfunc(spins[j], s_new) - bonds[n].bondfunc(spins[j], spins[m]);
                    }

                    if (randf() < 1. - std::exp(-dE/temperature)) {
                        c.push(j);
                    }
                }
            }

            if (is_ghost) {
                s0 = s0_new;
            } else {
                set_spin(m, s_new);
            }
        }
    }
}

void Spin3DModel::generate_mutation() {
    cluster_update();
}
#else
void Spin3DModel::metropolis_mutation() {
    if (acceptance > 0.5) {
        sigma = std::min(2., 1.01*sigma);

    } else {
        sigma = std::max(0.05, 0.99*sigma);
    }

    // Randomly generate mutation
    Eigen::Vector3d Gamma;
    Gamma << dist->sample(), dist->sample(), dist->sample();
    Eigen::Vector3d S2 = (spins[mut.i] + this->sigma*Gamma).normalized();


    // Store mutation for consideration
    this->mut.dS = S2 - spins[mut.i];
}

void Spin3DModel::generate_mutation() {
    mut.i = rand() % V;
    metropolis_mutation();
}
#endif


void Spin3DModel::accept_mutation() {
    return;
}

void Spin3DModel::reject_mutation() {
#ifndef CLUSTER_UPDATE
    set_spin(mut.i, spins[mut.i] - mut.dS);
#endif
}

double Spin3DModel::onsite_energy(int i) const {
#ifdef CLUSTER_UPDATE
    return onsite_func(s0.transpose()*spins[i]);
#else
    return onsite_func(spins[i]);
#endif
}

double Spin3DModel::bond_energy(int i) const {
    double E = 0.;
    int j;
    for (int n = 0; n < bonds.size(); n++) {
        j = neighbors[i][n];
        E += 0.5*bonds[n].bondfunc(spins[i], spins[j]);
    }

    return E;
}

double Spin3DModel::energy() const {
    double E = 0;

    for (int i = 0; i < V; i++) {
        E += onsite_energy(i);
        E += bond_energy(i);
    }

    return E;
}

double Spin3DModel::energy_change() {
#ifdef CLUSTER_UPDATE
    return -1.;
#else
    double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
    set_spin(mut.i, spins[mut.i] + mut.dS);
    double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

    return E2 - E1;
#endif
}

// Saves current spin configuration
void Spin3DModel::save_spins(std::string filename) {
    std::ofstream output_file;
    output_file.open(filename);
    output_file << N1 << "\t" << N2 << "\t" << N3 << "\t" << sl << "\n";
    int i; Eigen::Vector3d S;
    for (int i = 0; i < V; i++) {
#ifdef CLUSTER_UPDATE
        S = s0.transpose()*spins[i];
#else
        S = spins[i];
#endif
        output_file << S[0] << "\t" << S[1] << "\t" << S[2];
        if (i < V-1) { output_file << "\t"; }
    }
    output_file.close();
}

std::map<std::string, Sample> Spin3DModel::take_samples() {
    std::map<std::string, Sample> samples;
    samples.emplace("energy", energy());
    samples.emplace("magnetization", get_magnetization().norm());
    return samples;
}
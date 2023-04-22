#include "Spin3DModel.h"
#include <stack>
#include <iostream>
#include <fstream>

#define DEFAULT_CLUSTER_UPDATE true

#define DEFAULT_SAMPLE_ENERGY true
#define DEFAULT_SAMPLE_MAGNETIZATION true
#define DEFAULT_SAMPLE_HELICITY false

#define DEFAULT_BOUNDARY_CONDITION "periodic"


#define GHOST -1

GaussianDist::GaussianDist(float mean, float std) {
	this->rd.seed(rand());
	this->gen = std::default_random_engine(rd());
	this->dist = std::normal_distribution<>(mean, std);
}

float GaussianDist::sample() {
	return dist(gen);
}

BoundaryCondition Spin3DModel::parse_boundary_condition(std::string s) {
    if (s == "periodic") return BoundaryCondition::Periodic;
    else if (s == "open") return BoundaryCondition::Open;
    else std::cout << "Invalid boundary condition: " << s << "\n";
    
    assert(false);
}

Spin3DModel::Spin3DModel(Params &params) : nsteps(0), accepted(0) {
    cluster_update = params.get<int>("cluster_update", DEFAULT_CLUSTER_UPDATE);
    
    sample_helicity = params.get<int>("sample_helicity", DEFAULT_SAMPLE_HELICITY);
    sample_magnetization = params.get<int>("sample_magnetization", DEFAULT_SAMPLE_MAGNETIZATION);
    sample_energy = params.get<int>("sample_energy", DEFAULT_SAMPLE_ENERGY);

    bcx = parse_boundary_condition(params.get<std::string>("bcx", DEFAULT_BOUNDARY_CONDITION));
    bcy = parse_boundary_condition(params.get<std::string>("bcy", DEFAULT_BOUNDARY_CONDITION));
    bcz = parse_boundary_condition(params.get<std::string>("bcz", DEFAULT_BOUNDARY_CONDITION));
}

void Spin3DModel::init_params(int sl, int N1, int N2=-1, int N3=-1) {
    this->sl = sl;
    this->N1 = N1;
    if (N2 == -1) { this->N2 = N1; } else { this->N2 = N2; }
    if (N3 == -1) { this->N3 = N1; } else { this->N3 = N3; }
    this->V = N1*N2*N3*sl;

    this->spins = std::vector<Eigen::Vector3d>(V);
    if (cluster_update)
        this->neighbors = std::vector<std::vector<Bond>>(V+1, std::vector<Bond>(0));
    else
        this->neighbors = std::vector<std::vector<Bond>>(V, std::vector<Bond>(0));
}

void Spin3DModel::init() {
    if (cluster_update) {
        // Connect every site to the ghost 
        for (int i = 0; i < V; i++) {
            neighbors[V].push_back(Bond{i, GHOST});
            neighbors[i].push_back(Bond{V, GHOST});
        }

        this->s0 = Eigen::Matrix3d::Identity();
    }

    this->randomize_spins();

    this->acceptance = 0.5;
    this->sigma = 0.25;

    this->dist = GaussianDist(0., 1.0);

    this->mut.i = 0;
}

void Spin3DModel::randomize_spins() {
    for (int i = 0; i < V; i++)
        spins[i] = Eigen::Vector3d::Random(3).normalized();
}

void Spin3DModel::add_bond(int d1, int d2, int d3, int ds, Eigen::Vector3d v, std::function<double(Eigen::Vector3d, Eigen::Vector3d)> bondfunc) {
    HeisBond b{d1, d2, d3, ds, v, bondfunc};
    this->bonds.push_back(b);
    int i; int j;

    for (uint i = 0; i < V; i++) {
        Eigen::Vector4i idx = tensor_idx(i);

        uint nx = idx[0] + b.d1;
        uint ny = idx[1] + b.d2;
        uint nz = idx[2] + b.d3;
        uint ns = idx[3] + b.ds;

        if (bcx == BoundaryCondition::Open) {
            if (nx < 0 || nx > N1) continue;
        } else if (bcx == BoundaryCondition::Periodic) {
            nx = mod(nx, N1);
        }

        if (bcy == BoundaryCondition::Open) {
            if (ny < 0 || ny > N2) continue;
        } else if (bcx == BoundaryCondition::Periodic) {
            ny = mod(ny, N2);
        }

        if (bcz == BoundaryCondition::Open) {
            if (nz < 0 || nz > N3) continue;
        } else if (bcz == BoundaryCondition::Periodic) {
            nz = mod(nz, N3);
        }

        uint j = flat_idx(nx, ny, nz, ns);
        neighbors[i].push_back(Bond{j, bonds.size() - 1});
    }

//    for (int n1 = 0; n1 < N1; n1++) {
//        for (int n2 = 0; n2 < N2; n2++) {
//            for (int n3 = 0; n3 < N3; n3++) {
//                for (int s = 0; s < sl; s++) {
//                    i = flat_idx(n1, n2, n3, s);
//                    j = flat_idx(mod(n1 + b.d1, N1), mod(n2 + b.d2, N2), mod(n3 + b.d3, N3), mod(s + b.ds, sl));
//                    neighbors[i].push_back(j);
//                }
//            }
//        }
//    }

    double f = v[0]*alpha;
    Eigen::Matrix3d R;
    R << std::cos(f), -std::sin(f), 0,
         std::sin(f),  std::cos(f), 0.,
         0.,           0.,          1.;
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
    Eigen::Vector3d S1 = get_spin(i);
    Eigen::Vector3d S2;
    for (auto const &[j, b] : neighbors[i]) {
        if (b == GHOST) continue;
        
        S2 = get_spin(j);

        E0 += bonds[b].bondfunc(S1, S2);

        E1 += bonds[b].bondfunc(S1, R1s[b]*S2);
        Em1 += bonds[b].bondfunc(S1, R1s[b].transpose()*S2);

        E2 += bonds[b].bondfunc(S1, R2s[b]*S2);
        Em2 += bonds[b].bondfunc(S1, R2s[b].transpose()*S2);

        E3 += bonds[b].bondfunc(S1, R3s[b]*S2);
        Em3 += bonds[b].bondfunc(S1, R3s[b].transpose()*S2);
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
        for (int j = 0; j < 4; j++)
            twist[j] += twist_i[j];
    }
    
    return twist;
}

Eigen::Vector3d Spin3DModel::get_magnetization() const {
    Eigen::Vector3d M = Eigen::Vector3d::Constant(0);
    for (int i = 0; i < V; i++)
        M += get_spin(i);
    
    return M/V;
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
    Eigen::Vector3d dSdX; dSdX << 0., 0., 0.;
    Eigen::Vector3d dSdY; dSdY << 0., 0., 0.;
    for (auto const &[j, b] : neighbors[i]) {
        if (b == GHOST) continue;

        if (bonds[b].v[0] != 0.)
            dSdX += bonds[b].v[0]*(spins[j] - spins[i]);

        if (bonds[b].v[1] != 0.)
            dSdX += bonds[b].v[1]*(spins[j] - spins[i]);

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

void Spin3DModel::cluster_mutation() {
    s.clear();

    std::stack<int> c;
    int m = rand() % V;
    c.push(m);

    Eigen::Vector3d ax; ax << dist.sample(), dist.sample(), dist.sample();
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() - 2*ax*ax.transpose()/std::pow(ax.norm(), 2);

    double dE;
    Eigen::Matrix3d s0_new;
    Eigen::Vector3d s_new;
    bool is_ghost; bool neighbor_is_ghost;
    while (!c.empty()) {
        m = c.top();
        c.pop();


        if (!s.count(m)) {
            s.insert(m);

            is_ghost = (m == V);
            if (is_ghost) // Site is ghost
                s0_new = R*s0;
            else
                s_new = R*spins[m];

            for (auto const &[j, b] : neighbors[m]) {
                neighbor_is_ghost = (j == V);

                if (!s.count(j)) {
                    if (neighbor_is_ghost)
                        dE = onsite_func(s0.transpose()*s_new) - onsite_func(s0.transpose()*spins[m]);
                    else if (is_ghost)
                        dE = onsite_func(s0_new.transpose()*spins[j]) - onsite_func(s0.transpose()*spins[j]);
                    else // Normal bond
                        dE = bonds[b].bondfunc(spins[j], s_new) - bonds[b].bondfunc(spins[j], spins[m]);

                    if (randf() < 1. - std::exp(-dE/temperature))
                        c.push(j);
                }
            }

            if (is_ghost)
                s0 = s0_new;
            else
                set_spin(m, s_new);
        }
    }
}

void Spin3DModel::metropolis_mutation() {
    nsteps++;
    acceptance = accepted/nsteps;
    if (acceptance > 0.5)
        sigma = std::min(2., 1.01*sigma);
    else
        sigma = std::max(0.05, 0.99*sigma);

    // Randomly generate mutation
    Eigen::Vector3d Gamma;
    Gamma << dist.sample(), dist.sample(), dist.sample();
    Eigen::Vector3d S2 = (spins[mut.i] + this->sigma*Gamma).normalized();


    // Store mutation for consideration
    this->mut.dS = S2 - spins[mut.i];
}

void Spin3DModel::generate_mutation() {
    if (cluster_update)
        cluster_mutation();
    else {
        mut.i = rand() % V;
        metropolis_mutation();
    }
}


void Spin3DModel::accept_mutation() {
    accepted++;
    return;
}

void Spin3DModel::reject_mutation() {
    if (!cluster_update)
        set_spin(mut.i, spins[mut.i] - mut.dS);
}

double Spin3DModel::onsite_energy(int i) const {
    if (cluster_update)
        return onsite_func(s0.transpose()*spins[i]);
    else
        return onsite_func(spins[i]);
}

double Spin3DModel::bond_energy(int i) const {
    double E = 0.;
    double dE;
    for (auto const &[j, b] : neighbors[i]) {
        if (b == GHOST) continue;
        dE = 0.5*bonds[b].bondfunc(spins[i], spins[j]);
        E += dE;
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
    if (cluster_update)
        return -1.;

    double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
    set_spin(mut.i, spins[mut.i] + mut.dS);
    double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

    return E2 - E1;
}

// Saves current spin configuration
void Spin3DModel::save_spins(std::string filename) {
    std::ofstream output_file;
    output_file.open(filename);
    output_file << N1 << "\t" << N2 << "\t" << N3 << "\t" << sl << "\n";
    int i; Eigen::Vector3d S;
    for (int i = 0; i < V; i++) {
        if (cluster_update)
            S = s0.transpose()*spins[i];
        else
            S = spins[i];

        output_file << S[0] << "\t" << S[1] << "\t" << S[2];
        if (i < V-1) { output_file << "\t"; }
    }
    output_file.close();
}

void Spin3DModel::add_helicity_samples(std::map<std::string, Sample> &samples) const {
    std::vector<double> tterms = twist_derivatives();
    samples.emplace("d1E", tterms[0]);
    samples.emplace("d2E", tterms[1]);
    samples.emplace("d3E", tterms[2]);
    samples.emplace("d4E", tterms[3]);
}

void Spin3DModel::add_magnetization_samples(std::map<std::string, Sample> &samples) const {
    Eigen::Vector3d m = get_magnetization();
    samples.emplace("mx", m[0]);
    samples.emplace("my", m[1]);
    samples.emplace("mz", m[2]);
    samples.emplace("magnetization", m.norm());
}

std::map<std::string, Sample> Spin3DModel::take_samples() {
    std::map<std::string, Sample> samples;

    if (sample_energy)
        samples.emplace("energy", energy());
    if (sample_magnetization)
        add_magnetization_samples(samples);
    if (sample_helicity)
        add_helicity_samples(samples);

    return samples;
}
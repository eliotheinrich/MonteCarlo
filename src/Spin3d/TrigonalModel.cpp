#include "TrigonalModel.h"

TrigonalModel::TrigonalModel(Params &params) {
    this->N = params.get<int>("system_size");
    this->L = params.get<int>("layers", DEFAULT_LAYERS);
    Spin3DModel::init_params(1, N, N, L);

    this->J1 = params.get<float>("J1");
    this->J2 = params.get<float>("J2");
    this->K1 = params.get<float>("K1");
    this->K2 = params.get<float>("K2");
    this->K3 = params.get<float>("K3");
    Eigen::Vector3d B; B << params.get<float>("Bx"), params.get<float>("By"), params.get<float>("Bz");
    this->B = B;

    this->R << std::sqrt(3)/2., -0.5, 0.,
                0.5, std::sqrt(3)/2., 0.,
                0., 0., 1.;
    this->mut_mode = 0;
    this->mut_counter = 0; 

    float J1 = this->J1;
    std::function<float(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc = 
    [J1](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
        return -J1*S1.dot(S2);
    };

    Eigen::Vector3d v1; v1 << 1., 0., 0.;
    Eigen::Vector3d v2; v2 << 0.5, std::sqrt(3)/2., 0.;
    Eigen::Vector3d v3; v3 << 0.5, -std::sqrt(3)/2., 0.;
    this->add_bond(1,0,0,0, v1, bondfunc);
    this->add_bond(-1,0,0,0, -v1, bondfunc);
    this->add_bond(0,1,0,0, v2, bondfunc);
    this->add_bond(0,-1,0,0, -v2, bondfunc);
    this->add_bond(1,-1,0,0, v3, bondfunc);
    this->add_bond(-1,1,0,0, -v3, bondfunc);

    float J2 = this->J2;
    std::function<float(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc_inter = 
    [J2](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
        return J2*S1.dot(S2);
    };

    if ((std::abs(J2) > 1e-6) && (L > 1)) {
        Eigen::Vector3d v4; v4 << 0., 0., 1.;
        this->add_bond(0,0,1,0, v4, bondfunc_inter);
        this->add_bond(0,0,-1,0, -v4, bondfunc_inter);
    }

    this->sample_layer_magnetization = params.get<int>("sample_layer_magnetization", DEFAULT_SAMPLE_LAYER_MAGNETIZATION);
    this->sample_intensity = params.get<int>("sample_intensity", DEFAULT_SAMPLE_INTENSITY);
}

Eigen::Vector3d TrigonalModel::molecular_field(int i) const {
    float x = spins[i][0];
    float y = spins[i][1];
    float z = spins[i][2];
    Eigen::Vector3d H; H << K3*(6*std::pow(x, 5) - 60*std::pow(x, 3)*std::pow(y, 2) + 30*x*std::pow(y, 4)),
                            K3*(-6*std::pow(y, 5) + 60*std::pow(x, 2)*std::pow(y, 3) - 30*std::pow(x, 4)*y),
                            K1*std::pow(spins[i][2], 2);
    H += B;

    int j;
    for (int n = 0; n < 6; n++) {
        j = neighbors[i][n];
        H -= J1*spins[j];
    }
    for (int n = 6; n < 8; n++) {
        j = neighbors[i][n];
        H += J2*spins[j];
    }

    return H;
}

void TrigonalModel::generate_mutation() {
#ifdef CLUSTER_UPDATE
    cluster_update(); 
#else
    mut.i = rand() % V;
    mut_counter++;

    if (mut_counter == V) {
        mut_counter = 0;
        mut_mode++;
    }

    if (mut_mode < 10) {
        over_relaxation_mutation();
    } else if (mut_mode < 14) {
        metropolis_mutation();
    } else {
        metropolis_mutation();
        mut_mode = 0;
    }
#endif
}

void TrigonalModel::over_relaxation_mutation() {
    Eigen::Vector3d H = B;

    int j;
    for (int n = 0; n < 6; n++) {
        j = neighbors[mut.i][n];
        H -= J1*spins[j];
    }

    if (bonds.size() > 7) {
        for (int n = 6; n < 8; n++) {
            j = neighbors[mut.i][n];
            H += J2*spins[j];
        }
    }


    this->mut.dS = -2*spins[mut.i] + 2.*spins[mut.i].dot(H)/std::pow(H.norm(),2) * H;
}

double TrigonalModel::onsite_func(const Eigen::Vector3d &S) const {
    double E = 0;

    // Onsite interactions
    E -= B.dot(S);

    float phi = std::atan2(S[1], S[0]);
    float theta;
    if (S[2] > 1.0) { theta = PI; }
    else if (S[2] < -1.0) { theta = -PI; }
    else { theta = std::acos(S[2]); }

    E += K1*S[2]*S[2];
    E += K2*std::pow(S[0]*S[0]+S[1]*S[1],2);
    E += K3*std::cos(6*phi)*std::pow(sin(theta), 6); // Sixfold magnetocrystalline field

    return E;
}

void TrigonalModel::rotate_spin(int i, Eigen::Vector3d v, float p) {
    float vx = v[0]/v.norm(); float vy = v[1]/v.norm(); float vz = v[2]/v.norm();
    Eigen::Matrix3d R; R << std::cos(p) + vx*vx*(1 - std::cos(p)), vx*vy*(1 - std::cos(p)) - vz*std::sin(p), vx*vz*(1 - std::cos(p)) + vy*std::sin(p),
                            vy*vx*(1 - std::cos(p)) + vz*std::sin(p), std::cos(p) + vy*vy*(1 - std::cos(p)), vy*vz*(1 - std::cos(p)) - vx*std::sin(p),
                            vz*vx*(1 - std::cos(p)) - vy*std::sin(p), vz*vy*(1 - std::cos(p)) + vx*std::sin(p), std::cos(p) + vz*vz*(1 - std::cos(p));
    spins[i] = R*spins[i];
}

void TrigonalModel::dynamic_step(float dt) {
    std::vector<Eigen::Vector3d> H = std::vector<Eigen::Vector3d>(V);
    int j;

    // Compute local molecular fields
    for (int i = 0; i < V; i++) {
        H[i] = molecular_field(i);
    }

    float dT;
    float Hm;

    // Precess around local molecular field
    for (int i = 0; i < V; i++) {
        Hm = H[i].norm();
        dT = Hm*dt;
        spins[i] = cos(dT)*spins[i] + sin(dT)*H[i].cross(spins[i])/Hm + (1 - cos(dT))*H[i].dot(spins[i])*H[i]/pow(Hm, 2);
    }
}

#define c_bond 7.7072
#define a_bond 4.6926

Eigen::Vector3d TrigonalModel::rel_pos(uint i) const {
    Eigen::Vector4i idx = tensor_idx(i);
    uint n1 = idx[0]; uint n2 = idx[1]; uint n3 = idx[2];

    Eigen::Vector3d pos; pos << a_bond*0.5*(n1 + n2), a_bond*std::sqrt(3)/2.*(n1 - n2), c_bond*n3;
    return pos;
}

double TrigonalModel::intensity(Eigen::Vector3d Q) const {
    Eigen::Vector3cd structure_factor; structure_factor << 0., 0., 0.;
    for (uint i = 0; i < V; i++) {
        structure_factor += spins[i]*std::exp(std::complex<double>(0., Q.dot(rel_pos(i))));
    }
    return std::pow(std::abs(structure_factor.norm())/V, 2);
}

void TrigonalModel::add_layer_magnetization_samples(std::map<std::string, Sample> &samples) const {
    std::vector<Eigen::Vector3d> magnetization(L, Eigen::Vector3d::Constant(0.));
    for (uint i = 0; i < V; i++) {
        auto idx = tensor_idx(i);
        uint layer = idx[2];
        magnetization[layer] += spins[i];
    }

    for (uint layer = 0; layer < L; layer++) {
        magnetization[layer] /= V/L;
        samples.emplace("magnetization_" + std::to_string(layer) + "x", magnetization[layer][0]);
        samples.emplace("magnetization_" + std::to_string(layer) + "y", magnetization[layer][1]);
        samples.emplace("magnetization_" + std::to_string(layer) + "z", magnetization[layer][2]);
    }
}

void TrigonalModel::add_intensity_samples(std::map<std::string, Sample> &samples) const {
    for (uint i = 0; i < 30; i++) {
        Eigen::Vector3d q; q << 0., 0., double(i)/30*2*PI/c_bond;
        samples.emplace("intensity_" + std::to_string(i), intensity(q));
    }
}

std::map<std::string, Sample> TrigonalModel::take_samples() {
    std::map<std::string, Sample> samples = Spin3DModel::take_samples();
    /*
    std::vector<double> tterms = twist_derivatives();
    samples.emplace("d1E", tterms[0]);
    samples.emplace("d2E", tterms[1]);
    samples.emplace("d3E", tterms[2]);
    samples.emplace("d4E", tterms[3]);
    */

    if (sample_intensity)
        add_intensity_samples(samples);

    if (sample_layer_magnetization)
        add_layer_magnetization_samples(samples);
    
    if (std::abs(J2) > 1e-6) {
        double adjacency = 0.;
        for (uint i = 0; i < V; i++) {
            uint j = neighbors[i][6]; // Upstairs neighbor
            adjacency += spins[i].dot(spins[j]);
        }
        samples.emplace("adjacency", adjacency/V);
    }


    return samples;
}


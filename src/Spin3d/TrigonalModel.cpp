#include "TrigonalModel.h"

#define DEFAULT_LAYERS 1

#define DEFAULT_SAMPLE_LAYER_MAGNETIZATION false

#define DEFAULT_SAMPLE_INTENSITYX false
#define DEFAULT_SAMPLE_INTENSITYY false
#define DEFAULT_SAMPLE_INTENSITYZ false
#define DEFAULT_MAX_L 1.
#define DEFAULT_MIN_L 0.
#define DEFAULT_INTENSITY_RESOLUTION 30

#define DEFAULT_FM_LAYERS 0


TrigonalModel::TrigonalModel(Params &params) : Spin3DModel(params) {
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

    this->fm_layers = params.get<int>("fm_layers", DEFAULT_FM_LAYERS);
    this->mut_mode = 0;
    this->mut_counter = 0; 

    std::function<float(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc = 
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
        return -this->J1*S1.dot(S2);
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

    std::function<float(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc_inter = 
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
        return this->J2*S1.dot(S2);
    };

    std::function<float(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc_inter_fm = 
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
        return -this->J2*S1.dot(S2);
    };


    std::function<bool(uint, uint)> fm_layer_filter =
    [this](uint i, uint j) {
        auto idx1 = tensor_idx(i);
        auto idx2 = tensor_idx(j);
        return idx1[2] < this->fm_layers && idx2[2] < this->fm_layers;
    };

    std::function<bool(uint, uint)> afm_layer_filter = [&fm_layer_filter](uint i, uint j) { return !fm_layer_filter(i, j); };

    Eigen::Vector3d v4; v4 << 0., 0., 1.;
    this->add_bond(0, 0, 1, 0, v4, bondfunc_inter,    afm_layer_filter);
    this->add_bond(0, 0,-1, 0,-v4, bondfunc_inter,    afm_layer_filter);
    this->add_bond(0, 0, 1, 0, v4, bondfunc_inter_fm, fm_layer_filter);
    this->add_bond(0, 0,-1, 0,-v4, bondfunc_inter_fm, fm_layer_filter);




    this->sample_layer_magnetization = params.get<int>("sample_layer_magnetization", DEFAULT_SAMPLE_LAYER_MAGNETIZATION);

    this->sample_intensityx = params.get<int>("sample_intensityx", DEFAULT_SAMPLE_INTENSITYX);
    this->sample_intensityy = params.get<int>("sample_intensityy", DEFAULT_SAMPLE_INTENSITYY);
    this->sample_intensityz = params.get<int>("sample_intensityz", DEFAULT_SAMPLE_INTENSITYZ);
    this->max_L = params.get<float>("max_L", DEFAULT_MAX_L);
    this->min_L = params.get<float>("min_L", DEFAULT_MIN_L);
    this->intensity_resolution = params.get<int>("intensity_resolution", DEFAULT_INTENSITY_RESOLUTION);
}

Eigen::Vector3d TrigonalModel::molecular_field(int i) const {
    auto spin = get_spin(i);
    float x = spin[0];
    float y = spin[1];
    float z = spin[2];
    Eigen::Vector3d H; H << K3*(6*std::pow(x, 5) - 60*std::pow(x, 3)*std::pow(y, 2) + 30*x*std::pow(y, 4)),
                            K3*(-6*std::pow(y, 5) + 60*std::pow(x, 2)*std::pow(y, 3) - 30*std::pow(x, 4)*y),
                            K1*std::pow(z, 2);
    H += B;

    for (uint n = 0; n < 6; n++) {
        auto [j, _] = neighbors[i][n];
        H -= J1*get_spin(j);
    }
    for (uint n = 6; n < 8; n++) {
        auto [j, _] = neighbors[i][n];
        H += J2*get_spin(j);
    }

    return H;
}

void TrigonalModel::generate_mutation() {
    if (cluster_update)
        cluster_mutation(); 
    else {
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
    }
}

void TrigonalModel::over_relaxation_mutation() {
    Eigen::Vector3d H = B;

    for (uint n = 0; n < 6; n++) {
        auto [j, _] = neighbors[mut.i][n];
        H -= J1*get_spin(j);
    }

    if (bonds.size() > 7) {
        for (uint n = 6; n < 8; n++) {
            auto [j, _] = neighbors[mut.i][n];
            H += J2*get_spin(j);
        }
    }


    this->mut.dS = -2*get_spin(mut.i) + 2.*get_spin(mut.i).dot(H)/std::pow(H.norm(),2) * H;
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
    set_spin(i, R*get_spin(i));
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
        auto spin = get_spin(i);
        set_spin(i, cos(dT)*spin + sin(dT)*H[i].cross(spin)/Hm + (1 - cos(dT))*H[i].dot(spin)*H[i]/pow(Hm, 2));
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
        structure_factor += get_spin(i)*std::exp(std::complex<double>(0., Q.dot(rel_pos(i))));
    }
    return std::pow(std::abs(structure_factor.norm())/V, 2);
}

void TrigonalModel::add_layer_magnetization_samples(std::map<std::string, Sample> &samples) const {
    std::vector<Eigen::Vector3d> magnetization(L, Eigen::Vector3d::Constant(0.));
    for (uint i = 0; i < V; i++) {
        auto idx = tensor_idx(i);
        uint layer = idx[2];
        magnetization[layer] += get_spin(i);
    }

    for (uint layer = 0; layer < L; layer++) {
        magnetization[layer] /= V/L;
        samples.emplace("magnetization_" + std::to_string(layer) + "x", magnetization[layer][0]);
        samples.emplace("magnetization_" + std::to_string(layer) + "y", magnetization[layer][1]);
        samples.emplace("magnetization_" + std::to_string(layer) + "z", magnetization[layer][2]);
    }
}

void TrigonalModel::add_intensityx_samples(std::map<std::string, Sample> &samples) const {
    for (uint i = 0; i < intensity_resolution; i++) {
        float L = (i*max_L + (intensity_resolution - i)*min_L)/intensity_resolution;
        float d = 2.*PI*L/a_bond;
        Eigen::Vector3d q; q << d, d/std::sqrt(3), 0;
        samples.emplace("intensityx_" + std::to_string(i), intensity(q));
    }
}

void TrigonalModel::add_intensityy_samples(std::map<std::string, Sample> &samples) const {
    for (uint i = 0; i < intensity_resolution; i++) {
        float L = (i*max_L + (intensity_resolution - i)*min_L)/intensity_resolution;
        Eigen::Vector3d q; q << 0, 4.*PI*L/(std::sqrt(3)*a_bond), 0;
        samples.emplace("intensityy_" + std::to_string(i), intensity(q));
    }
}

void TrigonalModel::add_intensityz_samples(std::map<std::string, Sample> &samples) const {
    for (uint i = 0; i < intensity_resolution; i++) {
        float L = (i*max_L + (intensity_resolution - i)*min_L)/intensity_resolution;
        Eigen::Vector3d q; q << 0., 0., 2.*PI*L/c_bond;
        samples.emplace("intensityz_" + std::to_string(i), intensity(q));
    }
}

std::map<std::string, Sample> TrigonalModel::take_samples() {
    std::map<std::string, Sample> samples = Spin3DModel::take_samples();
    
    if (sample_intensityx)
        add_intensityx_samples(samples);
    if (sample_intensityy)
        add_intensityy_samples(samples);
    if (sample_intensityz)
        add_intensityz_samples(samples);

    if (sample_layer_magnetization)
        add_layer_magnetization_samples(samples);
    
    return samples;
}


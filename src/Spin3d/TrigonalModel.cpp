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


TrigonalModel::TrigonalModel(dataframe::Params &params, uint32_t num_threads) : Spin3DModel(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");
  L = dataframe::utils::get<int>(params, "layers", DEFAULT_LAYERS);

  J1 = dataframe::utils::get<double>(params, "J1");
  J2 = dataframe::utils::get<double>(params, "J2");
  J3 = dataframe::utils::get<double>(params, "J3", J2);
  K1 = dataframe::utils::get<double>(params, "K1");
  K2 = dataframe::utils::get<double>(params, "K2");
  K3 = dataframe::utils::get<double>(params, "K3");
  Eigen::Vector3d Bt; Bt << dataframe::utils::get<double>(params, "Bx"), dataframe::utils::get<double>(params, "By"), dataframe::utils::get<double>(params, "Bz");
  B = Bt;

  fm_layers = dataframe::utils::get<int>(params, "fm_layers", DEFAULT_FM_LAYERS);

  sample_layer_magnetization = dataframe::utils::get<int>(params, "sample_layer_magnetization", DEFAULT_SAMPLE_LAYER_MAGNETIZATION);

  sample_intensityx = dataframe::utils::get<int>(params, "sample_intensityx", DEFAULT_SAMPLE_INTENSITYX);
  sample_intensityy = dataframe::utils::get<int>(params, "sample_intensityy", DEFAULT_SAMPLE_INTENSITYY);
  sample_intensityz = dataframe::utils::get<int>(params, "sample_intensityz", DEFAULT_SAMPLE_INTENSITYZ);
  max_L = dataframe::utils::get<double>(params, "max_L", DEFAULT_MAX_L);
  min_L = dataframe::utils::get<double>(params, "min_L", DEFAULT_MIN_L);
  intensity_resolution = dataframe::utils::get<int>(params, "intensity_resolution", DEFAULT_INTENSITY_RESOLUTION);

  mut_mode = 0;
  mut_counter = 0; 

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc = 
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return -J1*S1.dot(S2);
    };

  Eigen::Vector3d v1; v1 << 1., 0., 0.;
  Eigen::Vector3d v2; v2 << 0.5, std::sqrt(3)/2., 0.;
  Eigen::Vector3d v3; v3 << 0.5, -std::sqrt(3)/2., 0.;
  add_bond(1,0,0,0, v1, bondfunc);
  add_bond(-1,0,0,0, -v1, bondfunc);
  add_bond(0,1,0,0, v2, bondfunc);
  add_bond(0,-1,0,0, -v2, bondfunc);
  add_bond(1,-1,0,0, v3, bondfunc);
  add_bond(-1,1,0,0, -v3, bondfunc);

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc_inter = 
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return J2*S1.dot(S2);
    };

  std::function<double(const Eigen::Vector3d &, const Eigen::Vector3d &)> bondfunc_inter_fm = 
    [this](const Eigen::Vector3d &S1, const Eigen::Vector3d &S2) {
      return -J3*S1.dot(S2);
    };


  std::function<bool(uint32_t, uint32_t)> fm_layer_filter = [this](uint32_t i, uint32_t j) {
    auto idx1 = tensor_idx(i);
    auto idx2 = tensor_idx(j);
    return idx1[2] < fm_layers && idx2[2] < fm_layers;
  };

  std::function<bool(uint32_t, uint32_t)> afm_layer_filter = [fm_layer_filter](uint32_t i, uint32_t j) {
    return !fm_layer_filter(i, j); 
  };

  Eigen::Vector3d v4; v4 << 0., 0., 1.;
  add_bond(0, 0, 1, 0, v4, bondfunc_inter,    afm_layer_filter);
  add_bond(0, 0,-1, 0,-v4, bondfunc_inter,    afm_layer_filter);
  add_bond(0, 0, 1, 0, v4, bondfunc_inter_fm, fm_layer_filter);
  add_bond(0, 0,-1, 0,-v4, bondfunc_inter_fm, fm_layer_filter);
  
  Spin3DModel::init(1, N, N, L);
}

Eigen::Vector3d TrigonalModel::molecular_field(uint32_t i) const {
  auto spin = get_spin(i);
  double x = spin[0];
  double y = spin[1];
  double z = spin[2];
  Eigen::Vector3d H; H << K3*(6*std::pow(x, 5) - 60*std::pow(x, 3)*std::pow(y, 2) + 30*x*std::pow(y, 4)),
                          K3*(-6*std::pow(y, 5) + 60*std::pow(x, 2)*std::pow(y, 3) - 30*std::pow(x, 4)*y),
                          K1*std::pow(z, 2);
  H += B;

  for (uint32_t n = 0; n < 6; n++) {
    auto [j, _] = neighbors[i][n];
    H -= J1*get_spin(j);
  }
  for (uint32_t n = 6; n < 8; n++) {
    auto [j, _] = neighbors[i][n];
    H += J2*get_spin(j);
  }

  return H;
}

void TrigonalModel::generate_mutation() {
  if (cluster_update) {
    cluster_mutation(); 
  } else {
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

// This function has a high capacity for bugs; beware when changing things
void TrigonalModel::over_relaxation_mutation() {
  Eigen::Vector3d H = B;

  for (auto const &[j, n] : neighbors[mut.i]) {
    if (n < 6) {
      H -= J1*get_spin(j);
    } else if (n < 8) {
      H += J2*get_spin(j);
    } else if (n < 10) {
      H -= J3*get_spin(j);
    }
  }

  mut.dS = -2*get_spin(mut.i) + 2.*get_spin(mut.i).dot(H)/std::pow(H.norm(),2) * H;
}

double TrigonalModel::onsite_func(const Eigen::Vector3d &S) const {
  double E = 0;

  // Onsite interactions
  E -= B.dot(S);

  double phi = std::atan2(S[1], S[0]);
  double theta;
  if (S[2] > 1.0) { 
    theta = PI; 
  } else if (S[2] < -1.0) { 
    theta = -PI; 
  } else { 
    theta = std::acos(S[2]);
  }

  E += K1*S[2]*S[2];
  E += K2*std::pow(S[0]*S[0]+S[1]*S[1],2);
  E += K3*std::cos(6*phi)*std::pow(sin(theta), 6); // Sixfold magnetocrystalline field

  return E;
}

void TrigonalModel::rotate_spin(uint32_t i, Eigen::Vector3d v, double p) {
  double vx = v[0]/v.norm(); double vy = v[1]/v.norm(); double vz = v[2]/v.norm();
  Eigen::Matrix3d R; R << std::cos(p) + vx*vx*(1 - std::cos(p)), vx*vy*(1 - std::cos(p)) - vz*std::sin(p), vx*vz*(1 - std::cos(p)) + vy*std::sin(p),
                          vy*vx*(1 - std::cos(p)) + vz*std::sin(p), std::cos(p) + vy*vy*(1 - std::cos(p)), vy*vz*(1 - std::cos(p)) - vx*std::sin(p),
                          vz*vx*(1 - std::cos(p)) - vy*std::sin(p), vz*vy*(1 - std::cos(p)) + vx*std::sin(p), std::cos(p) + vz*vz*(1 - std::cos(p));
  set_spin(i, R*get_spin(i));
}

void TrigonalModel::dynamic_step(double dt) {
  std::vector<Eigen::Vector3d> H = std::vector<Eigen::Vector3d>(V);

  // Compute local molecular fields
  for (uint32_t i = 0; i < V; i++) {
    H[i] = molecular_field(i);
  }

  double dT;
  double Hm;

  // Precess around local molecular field
  for (uint32_t i = 0; i < V; i++) {
    Hm = H[i].norm();
    dT = Hm*dt;
    auto spin = get_spin(i);
    set_spin(i, cos(dT)*spin + sin(dT)*H[i].cross(spin)/Hm + (1 - cos(dT))*H[i].dot(spin)*H[i]/pow(Hm, 2));
  }
}

#define c_bond 7.7072
#define a_bond 4.6926

Eigen::Vector3d TrigonalModel::rel_pos(uint32_t i) const {
  Eigen::Vector4i idx = tensor_idx(i);
  uint32_t n1 = idx[0]; 
  uint32_t n2 = idx[1]; 
  uint32_t n3 = idx[2];

  Eigen::Vector3d pos; pos << a_bond*0.5*(n1 + n2), a_bond*std::sqrt(3)/2.*(n1 - n2), c_bond*n3;
  return pos;
}

double TrigonalModel::intensity(Eigen::Vector3d Q) const {
  Eigen::Vector3cd structure_factor; structure_factor << 0., 0., 0.;
  for (uint32_t i = 0; i < V; i++) {
    structure_factor += get_spin(i)*std::exp(std::complex<double>(0., Q.dot(rel_pos(i))));
  }

  return std::pow(std::abs(structure_factor.norm())/V, 2);
}

void TrigonalModel::add_layer_magnetization_samples(dataframe::data_t &samples) const {
  std::vector<Eigen::Vector3d> magnetization(L, Eigen::Vector3d::Constant(0.));
  for (uint32_t i = 0; i < V; i++) {
    auto idx = tensor_idx(i);
    uint layer = idx[2];
    magnetization[layer] += get_spin(i);
  }

  for (uint32_t layer = 0; layer < L; layer++) {
    magnetization[layer] /= V/L;
    samples.emplace("magnetization_" + std::to_string(layer) + "x", magnetization[layer][0]);
    samples.emplace("magnetization_" + std::to_string(layer) + "y", magnetization[layer][1]);
    samples.emplace("magnetization_" + std::to_string(layer) + "z", magnetization[layer][2]);
  }
}

void TrigonalModel::add_intensityx_samples(dataframe::data_t &samples) const {
  std::vector<double> intensity_samples(intensity_resolution);
  for (uint32_t i = 0; i < intensity_resolution; i++) {
    double L = (i*max_L + (intensity_resolution - i)*min_L)/intensity_resolution;
    double d = 2.*PI*L/a_bond;
    Eigen::Vector3d q; q << d, d/std::sqrt(3), 0;
    intensity_samples[i] = intensity(q);
  }

  samples.emplace("intensityy", intensity_samples);
}

void TrigonalModel::add_intensityy_samples(dataframe::data_t &samples) const {
  std::vector<double> intensity_samples(intensity_resolution);
  for (uint32_t i = 0; i < intensity_resolution; i++) {
    double L = (i*max_L + (intensity_resolution - i)*min_L)/intensity_resolution;
    Eigen::Vector3d q; q << 0, 4.*PI*L/(std::sqrt(3)*a_bond), 0;
    intensity_samples[i] = intensity(q);
  }

  samples.emplace("intensityy", intensity_samples);
}

void TrigonalModel::add_intensityz_samples(dataframe::data_t &samples) const {
  std::vector<double> intensity_samples(intensity_resolution);
  for (uint32_t i = 0; i < intensity_resolution; i++) {
    double L = (i*max_L + (intensity_resolution - i)*min_L)/intensity_resolution;
    Eigen::Vector3d q; q << 0., 0., 2.*PI*L/c_bond;
    intensity_samples[i] = intensity(q);
  }

  samples.emplace("intensityz", intensity_samples);
}

dataframe::data_t TrigonalModel::take_samples() {
  dataframe::data_t samples = Spin3DModel::take_samples();

  if (sample_intensityx) {
    add_intensityx_samples(samples);
  }

  if (sample_intensityy) {
    add_intensityy_samples(samples);
  }

  if (sample_intensityz) {
    add_intensityz_samples(samples);
  }

  if (sample_layer_magnetization) {
    add_layer_magnetization_samples(samples);
  }

  return samples;
}


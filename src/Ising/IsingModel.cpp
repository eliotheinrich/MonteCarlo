#include "IsingModel.h"
#include <fstream>

void IsingModel::init(uint32_t N1, uint32_t N2, uint32_t N3) {
  this->N1 = N1;
  this->N2 = N2;
  this->N3 = N3;
  V = N1*N2*N3;

  acceptance = 0.5;

  spins = std::vector<double>(V);
  randomize_spins();
}

void IsingModel::randomize_spins() {
  for (uint32_t i = 0; i < V; i++) {
    // For each site, initialize spin randomly
    if (rand() % 2) {
      spins[i] = 1.;
    } else {
      spins[i] = -1.;
    }
  }
}

double IsingModel::get_magnetization() const {
  double M = 0;
  for (uint32_t i = 0; i < V; i++) {
    M += spins[i];
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
  spins[mut.i] = -spins[mut.i];
}

double IsingModel::energy() const {
  double E = 0;

  for (uint32_t i = 0; i < V; i++) {
    E += onsite_energy(i);
    E += bond_energy(i);
  }
  return E;
}

double IsingModel::energy_change() {
  double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
  spins[mut.i] = -spins[mut.i];
  double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

  return E2 - E1;
}

uint64_t IsingModel::system_size() const {
  return V;
}

uint32_t IsingModel::flat_idx(uint32_t n1, uint32_t n2, uint32_t n3) const {
  return n1 + N1*(n2 + N2*n3);
}

Eigen::Vector3i IsingModel::tensor_idx(uint32_t i) const {
  uint32_t n1 = i % N1;
  i = i / N1;
  uint32_t n2 = i % N2;
  i = i / N2;
  uint32_t n3 = i % N3;
  Eigen::Vector3i v; v << n1, n2, n3;
  return v;
}

dataframe::data_t IsingModel::take_samples() {
  dataframe::data_t samples;

  double E = energy();
  double M = get_magnetization();

  samples.emplace("energy", E);
  samples.emplace("energy_second_moment", E*E);
  samples.emplace("energy_third_moment", E*E*E);
  samples.emplace("energy_fourth_moment", E*E*E*E);

  samples.emplace("magnetization", std::abs(M));
  samples.emplace("signed_magnetization", M);
  samples.emplace("magnetization_second_moment", M*M);
  samples.emplace("magnetization_third_moment", M*M*M);
  samples.emplace("magnetization_fourth_moment", M*M*M*M);


  return samples;
}

// Saves current spin configuration
void IsingModel::save_spins(const std::string& filename) {
  std::ofstream output_file;
  output_file.open(filename);
  output_file << N1 << "\t" << N2 << "\t" << N3 << "\n";
  for (uint32_t i = 0; i < V; i++) {
    output_file << spins[i];
    if (i < V-1) { 
      output_file << "\t"; 
    }
  }
  output_file.close();
}        

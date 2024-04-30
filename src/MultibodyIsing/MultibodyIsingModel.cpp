#include "MultibodyIsingModel.h"

void MultibodyIsingModel::init(size_t L) {
  this->L = L;
  V = L*L;

  acceptance = 0.5;

  spins = std::vector<int>(V);
  randomize_spins();
  terms = std::vector<std::vector<MultibodyIsingTerm>>(V);
  onsite_potential = std::vector<double>(V, 0.0);
}

void MultibodyIsingModel::add_term(const std::vector<size_t>& inds, double J) {
  if (inds.size() == 1) {
    onsite_potential[inds[0]] += J;
  } else {
    MultibodyIsingTerm term = {J, inds};
    for (auto j : inds) {
      terms[j].push_back(term);
    }
  }
}

void MultibodyIsingModel::randomize_spins() {
  for (uint32_t i = 0; i < V; i++) {
    // For each site, initialize spin randomly
    if (rand() % 2) {
      spins[i] = 1;
    } else {
      spins[i] = -1;
    }
  }
}

double MultibodyIsingModel::get_magnetization() const {
  double M = 0;
  for (uint32_t i = 0; i < V; i++) {
    M += static_cast<double>(spins[i]);
  }

  return M/V;
}

void MultibodyIsingModel::generate_mutation() {
  // Randomly select a site to mutate
  mut.i = rand() % V;
}

void MultibodyIsingModel::accept_mutation() {
  return;
}

void MultibodyIsingModel::reject_mutation() {
  spins[mut.i] = -spins[mut.i];
}

double MultibodyIsingModel::onsite_energy(uint32_t i) const {
  return onsite_potential[i];
}

double MultibodyIsingModel::bond_energy(uint32_t i) const {
  double E = 0.0;
  for (auto const& term : terms[i]) {
    double s = 1.0;
    for (auto const j : term.inds) {
      s = s * spins[j];
    }
    E += term.J * s;
  }

  return E/2.0;
}

double MultibodyIsingModel::energy() const {
  double E = 0;

  for (uint32_t i = 0; i < V; i++) {
    E += onsite_energy(i);
    E += bond_energy(i);
  }
  return E;
}

double MultibodyIsingModel::energy_change() {
  double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
  spins[mut.i] = -spins[mut.i];
  double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

  return E2 - E1;
}

uint64_t MultibodyIsingModel::system_size() const {
  return V;
}

size_t MultibodyIsingModel::idx(size_t x, size_t y) const {
  return mod(x, L) + L * mod(y, L);
}

std::pair<int, int> MultibodyIsingModel::coordinates(size_t i) const {
  int x = i / L;
  int y = i % L;
  return std::make_pair(x, y);
}

dataframe::data_t MultibodyIsingModel::take_samples() {
  dataframe::data_t samples;

  double E = energy();
  double M = get_magnetization();

  samples.emplace("energy", E);

  samples.emplace("magnetization", std::abs(M));

  return samples;
}

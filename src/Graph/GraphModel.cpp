#include "GraphModel.h"

GraphModel::GraphModel(dataframe::ExperimentParams &params, uint32_t num_threads) : MonteCarloSimulator(params, num_threads) {
  edges = std::vector<std::vector<uint32_t>>(N, std::vector<uint32_t>(0));
  vals = std::vector<int>(N);

  for (uint32_t i = 0; i < N; i++) {
    vals[i] = 0;
    for (uint32_t j = 0; j < N; j++) {
      edges[i][j] = 0;
    }
  }
}

void GraphModel::init(uint64_t N) {
  this->N = N;
}

void GraphModel::generate_mutation() {
  mut.i = rand() % N;
  mut.j = rand() % N;
  while (mut.j == mut.i) {
    mut.j = rand() % N;
  }
}

void GraphModel::accept_mutation() {
  return;
}

void GraphModel::reject_mutation() {
  toggle_edge(mut.i, mut.j);
}

double GraphModel::energy() const {
  double E = 0;

  for (uint32_t i = 0; i < N; i++) {
    E += onsite_energy(i);
    E += bond_energy(i);
  }
  return E;
}

double GraphModel::energy_change() {
  double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
  toggle_edge(mut.i, mut.j);
  double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

  return E2 - E1;
}

void GraphModel::toggle_edge(uint32_t i, uint32_t j) {
  edges[i][j] = 1 - edges[i][j];
  edges[j][i] = 1 - edges[j][i];
}

uint32_t GraphModel::deg(uint32_t i) const {
  uint32_t d = 0;
  for (uint32_t j = 0; j < N; j++) {
    d += edges[i][j];
  }
  return d;
}

double GraphModel::get_connectivity() const {
  double c = 0.;
  for (uint32_t i = 0; i < N; i++) {
    c += deg(i);
  }
  return c/N;
}

dataframe::SampleMap GraphModel::take_samples() const {
  dataframe::SampleMap samples;

  dataframe::utils::emplace(samples, "connectivity", get_connectivity());

  return samples;
}

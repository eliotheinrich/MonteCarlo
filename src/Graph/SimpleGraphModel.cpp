#include "SimpleGraphModel.h"

SimpleGraphModel::SimpleGraphModel(dataframe::ExperimentParams &params, uint32_t num_threads) : GraphModel(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");
  J = dataframe::utils::get<double>(params, "J");
  GraphModel::init(N);
}

double SimpleGraphModel::onsite_energy(uint32_t i) const {
  double E = 0;

  return E;
}

double SimpleGraphModel::bond_energy(uint32_t i) const {
  double E = 0;
  E += deg(i);
  return J*E;
}

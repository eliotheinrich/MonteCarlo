#include "SquareIsingModel.h"

SquareIsingModel::SquareIsingModel(dataframe::Params &params, uint32_t num_threads) : IsingModel(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");
  L = dataframe::utils::get<int>(params, "layers", DEFAULT_LAYERS);

  J = dataframe::utils::get<double>(params, "J");
  B = dataframe::utils::get<double>(params, "B");

  IsingModel::init(N, N, L);
}

double SquareIsingModel::onsite_energy(uint32_t i) const {
  // Onsite interactions
  return -B*spins[i];
}

double SquareIsingModel::bond_energy(uint32_t i) const {
  double E = 0;

  Eigen::Vector3i idxs = tensor_idx(i);
  int n1 = idxs[0]; 
  int n2 = idxs[1]; 
  int n3 = idxs[2]; 

  // NN interactions
  E -= J*spins[i]*spins[flat_idx(mod(n1+1, N), n2, n3)];
  E -= J*spins[i]*spins[flat_idx(mod(n1-1, N), n2, n3)];

  E -= J*spins[i]*spins[flat_idx(n1, mod(n2+1, N), n3)];
  E -= J*spins[i]*spins[flat_idx(n1, mod(n2-1, N), n3)];

  return 0.5*E;
}

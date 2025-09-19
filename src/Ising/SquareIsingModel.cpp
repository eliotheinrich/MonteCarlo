#include "SquareIsingModel.h"

SquareIsingModel::SquareIsingModel(Params& params, uint32_t num_threads) : IsingModel(params, num_threads) {
  N = get<int>(params, "system_size");
  L = get<int>(params, "layers", DEFAULT_LAYERS);

  J = get<double>(params, "J");
  B = get<double>(params, "B");

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

Texture SquareIsingModel::get_texture() const {
  Texture texture(N, N);
  for (size_t i = 0; i < V; i++) {
    auto idxs = tensor_idx(i);
    size_t x = idxs[0];
    size_t y = idxs[1];
    if (spins[i] == 1) {
      texture.set(x, y, {1.0, 1.0, 1.0, 1.0});
    } else {
      texture.set(x, y, {0.0, 0.0, 0.0, 1.0});
    }
  }

  return texture;
}

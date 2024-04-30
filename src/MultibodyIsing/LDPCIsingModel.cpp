#include "LDPCIsingModel.h"

LDPCIsingModel::LDPCIsingModel(dataframe::Params &params, uint32_t num_threads) : MultibodyIsingModel(params, num_threads) {
  L = dataframe::utils::get<int>(params, "system_size");
  J = dataframe::utils::get<double>(params, "J", 1.0);
  impurity = dataframe::utils::get<double>(params, "impurity", 0.0);

  MultibodyIsingModel::init(L);

  for (size_t i = 0; i < V; i++) {
    auto [x, y] = coordinates(i);
    
    size_t i1, i2;
    if (x % 2 == y % 2) {
      i1 = idx(x+1, y);
      i2 = idx(x-1, y);
    } else {
      i1 = idx(x, y+1);
      i2 = idx(x, y-1);
    }

    std::vector<size_t> inds{i};
    if (randf() < 1.0 - impurity) {
      inds.push_back(i1);
      inds.push_back(i2);
    }

    add_term(inds, J);
  }
}

double LDPCIsingModel::onsite_energy(uint32_t i) const {
  return 0.0;
}

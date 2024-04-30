#include "SquareClockModel.h"

template <uint32_t q>
SquareClockModel<q>::SquareClockModel(dataframe::Params &params, uint32_t num_threads) : ClockModel<q>(params, num_threads) {
  N = dataframe::utils::get<int>(params, "system_size");
  L = dataframe::utils::get<int>(params, "layers", DEFAULT_LAYERS);

  J = dataframe::utils::get<double>(params, "J");

  for (uint32_t i = 0; i < q; i++) {
    for (uint32_t j = 0; j < q; j++) {
      bond_table[i][j] = -J*cos(2*PI/q*(i - j));
    }
  }

  std::function<double(uint32_t, uint32_t)> bondfunc = [this](uint32_t p1, uint32_t p2) {
    return this->bond_table[p1][p2];
  };



  Eigen::Vector3d v1; v1 << 1.,0.,0.;
  Eigen::Vector3d v2; v2 << 0.,1.,0.;
  this->add_bond(1,0,0,   v1, bondfunc);
  this->add_bond(-1,0,0, -v1, bondfunc);
  this->add_bond(0,1,0,   v2, bondfunc);
  this->add_bond(0,-1,0, -v2, bondfunc);

  ClockModel<q>::init(N, N, L);
}

template <uint32_t q>
double SquareClockModel<q>::onsite_energy(uint32_t i) const {
  return 0.;
}


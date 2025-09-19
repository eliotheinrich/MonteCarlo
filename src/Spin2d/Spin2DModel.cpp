#include "Spin2DModel.h"

#include <math.h>

Spin2DModel::Spin2DModel(Params &params, uint32_t num_threads) : MonteCarloSimulator(params, num_threads) {
  cluster_update = get<int>(params, "cluster_update", true);
}

void Spin2DModel::init(const Lattice<Spin2D>& lattice) {
  std::srand(randi());

  this->lattice = lattice;
  V = lattice.system_size();

  s0 = Eigen::Matrix2d::Identity();

  randomize_spins();

  mut.i = 0;
  sigma = 0.25;
}

void Spin2DModel::randomize_spins() {
  for (uint32_t i = 0; i < V; i++) {
    double p = 2*M_PI*randf();
    Spin2D s; s << std::cos(p), std::sin(p);
    set_spin(i, s);
  }
}

std::vector<double> Spin2DModel::twist_stiffness() const {
  // Returns the first and second derivative in response to a phase twist
  double E0 = 0.;
  double E1 = 0.;
  double E2 = 0.;
  double E3 = 0.;
  double Em1 = 0.;
  double Em2 = 0.;
  double Em3 = 0.;

  auto [R1s, R2s, R3s] = lattice.get_twist_matrices(alpha);

  for (uint32_t i = 0; i < V; i++) {
    for (auto const &[j, b] : lattice.neighbors[i]) {
      Spin2D S1 = get_spin(i);
      Spin2D S2 = get_spin(j);

      E0 += lattice.bonds[b].bondfunc(S1, S2);

      E1  += lattice.bonds[b].bondfunc(S1, R1s[b]*S2);
      Em1 += lattice.bonds[b].bondfunc(S1, R1s[b].transpose()*S2);

      E2  += lattice.bonds[b].bondfunc(S1, R2s[b]*S2);
      Em2 += lattice.bonds[b].bondfunc(S1, R2s[b].transpose()*S2);

      E3  += lattice.bonds[b].bondfunc(S1, R3s[b]*S2);
      Em3 += lattice.bonds[b].bondfunc(S1, R3s[b].transpose()*S2);
    }
  }

  // Compute derivates from finite difference
  double d1E = (1./12.*Em2 - 2./3.*Em1 + 2./3.*E1 - 1./12.*E2)/alpha/2.;
  double d2E = (-1./12.*Em2 + 4./3.*Em1 - 5./2.*E0 + 4./3.*E1 - 1./12.*E2)/pow(alpha, 2)/2.;
  double d3E = (1./8.*Em3 - 1.*Em2 + 13./8.*Em1 - 13./8.*E1 + 1.*E2 - 1./8.*E3)/pow(alpha, 3)/2.;
  double d4E = (-1./6.*Em3 + 2.*Em2 - 13./2.*Em1 + 28./3.*E0 - 13./2.*E1 + 2.*E2 - 1./6.*E3)/pow(alpha, 4)/2.;

  return std::vector<double>{d4E, d3E, d1E, pow(d2E,2), d2E, d1E*d3E, pow(d1E,2)*d2E, 
    d1E*d2E, pow(d1E,2), pow(d1E,4), pow(d1E,3)};
}

Eigen::Vector2d Spin2DModel::get_magnetization() const {
  Eigen::Vector2d M = Eigen::Vector2d::Zero();
  for (uint32_t i = 0; i < V; i++) {
    M += lattice.spins[i];
  }

  if (cluster_update) {
    return s0.transpose()*M/V;
  } else {
    return M/V;
  }
}

void Spin2DModel::generate_mutation() {
  if (cluster_update) {
    cluster_mutation();
  } else {
    mut.i = rand() % V;
    metropolis_mutation();
  }
}

void Spin2DModel::cluster_mutation() {
  s.clear();

  std::stack<uint32_t> c;
  uint32_t m = rand() % V;
  Eigen::Matrix2d s0_new;
  c.push(m);

  double p = (double) 2*M_PI*rand()/RAND_MAX;
  Eigen::Vector2d ax; ax << std::cos(p), std::sin(p);
  Eigen::Matrix2d R = Eigen::Matrix2d::Identity() - 2*ax*ax.transpose()/std::pow(ax.norm(), 2);

  Eigen::Vector2d s_new;
  while (!c.empty()) {
    m = c.top();
    c.pop();

    if (!s.count(m)) {
      s.insert(m);
      bool is_ghost = (m == V);
      if (is_ghost) { // Site is ghost
        s0_new = R*s0;
      } else {
        s_new = R*lattice.spins[m];
      }

      for (auto const &[j, b] : lattice.neighbors[m]) {
        if (!s.count(j)) {
          bool neighbor_is_ghost = (j == V);
          
          double dE;
          if (neighbor_is_ghost) {
            dE = onsite_func(s0.inverse()*s_new) - onsite_func(s0.inverse()*lattice.spins[m]);
          } else if (is_ghost) {
            dE = onsite_func(s0_new.inverse()*lattice.spins[j]) - onsite_func(s0.inverse()*lattice.spins[j]);
          } else { // Normal bond
            dE = lattice.bonds[b].bondfunc(lattice.spins[j], s_new) - lattice.bonds[b].bondfunc(lattice.spins[j], lattice.spins[m]);
          }

          if (randf() < 1. - std::exp(-dE/temperature)) {
            c.push(j);
          }
        }
      }

      if (is_ghost) {
        s0 = s0_new;
      } else {
        lattice.spins[m] = s_new;
      }
    }
  }
}

void Spin2DModel::metropolis_mutation() {
  double dp = sigma*(randf() - 0.5)*2.*M_PI;
  Eigen::Vector2d S1 = lattice.spins[mut.i];
  Eigen::Vector2d S2; S2 << std::cos(dp)*S1[0] - std::sin(dp)*S1[1],
                            std::cos(dp)*S1[1] + std::sin(dp)*S1[0];

  // Store mutation for consideration
  mut.dS = S2 - S1;
}


void Spin2DModel::accept_mutation() {
  return;
}

void Spin2DModel::reject_mutation() {
  if (!cluster_update) {
    lattice.spins[mut.i] = lattice.spins[mut.i] - mut.dS;
  }
}

double Spin2DModel::onsite_energy(int i) const {
  if (cluster_update) {
    return onsite_func(s0.transpose()*lattice.spins[i]);
  } else {
    return onsite_func(lattice.spins[i]);
  }
}

double Spin2DModel::bond_energy(int i) const {
  double E = 0.;
  for (auto const &[j, b] : lattice.neighbors[i]) {
    if (b == GHOST) {
      continue;
    }

    E += 0.5*lattice.bonds[b].bondfunc(lattice.spins[i], lattice.spins[j]);
  }

  return E;
}

double Spin2DModel::energy() const {
  double E = 0;

  for (uint32_t i = 0; i < V; i++) {
    E += onsite_energy(i);
    E += bond_energy(i);
  }

  return E;
}

double Spin2DModel::energy_change() {
  if (cluster_update) {
    return -1.;
  }

  double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
  lattice.spins[mut.i] = lattice.spins[mut.i] + mut.dS;
  double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

  return E2 - E1;
}

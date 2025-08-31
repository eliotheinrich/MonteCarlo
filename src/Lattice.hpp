#pragma once

#include <stdint.h>
#include <functional>
#include <optional>
#include <vector>

#include <Eigen/Dense>

#include <Graph.hpp>

constexpr uint32_t GHOST = -1;

struct LatticeDimension {
  size_t N;
  BoundaryCondition bc;
  Eigen::Vector3d v;
};

using SiteFilter = std::function<bool(uint32_t, uint32_t, uint32_t, uint32_t)>;

template <typename Spin>
using BondFunc = std::function<double(const Spin&, const Spin&)>;

using Spin2D = Eigen::Vector2d;
using Spin3D = Eigen::Vector3d;

template <typename Spin>
struct SpinBond {
  int d1;
  int d2;
  int d3;
  int ds;
  BondFunc<Spin> bondfunc;
  SiteFilter filter;

  SpinBond(int d1, int d2, int d3, int ds, BondFunc<Spin> bondfunc, std::optional<SiteFilter> filter=std::nullopt) : d1(d1), d2(d2), d3(d3), ds(ds), bondfunc(bondfunc) {
    if (filter) {
      this->filter = filter.value();
    } else {
      SiteFilter include_all = [](uint32_t, uint32_t, uint32_t, uint32_t) { return true; };
      this->filter = include_all;
    }
  }
};

using LatticeGraph = DirectedGraph<std::vector<double>, int>;

template <typename Spin>
struct Lattice {
  LatticeDimension dx;
  LatticeDimension dy;
  LatticeDimension dz;
  size_t num_sublattices;

  std::vector<Eigen::Vector3d> sublattice_vectors;

  std::vector<Spin> spins;
  std::vector<SpinBond<Spin>> bonds;
  std::vector<std::vector<Bond>> neighbors;

  size_t system_size() const {
    return dx.N * dy.N * dz.N * num_sublattices;
  }

  LatticeGraph to_graph() const {
    size_t num_spins = system_size();
    LatticeGraph graph(num_spins);

    for (size_t i = 0; i < num_spins; i++) {
      auto pos = position(i);
      std::vector<double> pos_vec(pos.begin(), pos.end());
      graph.set_val(i, pos_vec);
      for (const auto& [j, n] : neighbors[i]) {
        if (j == GHOST || j >= i) {
          continue;
        }

        graph.add_edge(i, j, n);
      }
    }

    return graph;
  }

  Lattice()=default;

  void add_ghost() {
    neighbors.push_back(std::vector<Bond>(0));

    // Connect every site to the ghost 
    size_t num_spins = spins.size();
    for (uint32_t i = 0; i < num_spins; i++) {
      neighbors[num_spins].push_back(Bond{i, GHOST});
      neighbors[i].push_back(Bond{num_spins, GHOST});
    }
  }

  Lattice(const LatticeDimension& dx, const LatticeDimension& dy, const LatticeDimension& dz,
    const std::vector<SpinBond<Spin>>& bonds,
    std::optional<std::vector<Eigen::Vector3d>> sublattice_vectors_opt=std::nullopt, 
    bool include_ghost_spin=false
  ) {
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;

    this->bonds = bonds;

    num_sublattices = 1;
    if (sublattice_vectors_opt) {
      this->sublattice_vectors = sublattice_vectors_opt.value();
      num_sublattices += this->sublattice_vectors.size();
    } else {
      this->sublattice_vectors = {};
    }

    size_t num_spins = system_size();
    spins = std::vector<Spin>(num_spins);
    neighbors = std::vector<std::vector<Bond>>(num_spins,   std::vector<Bond>(0));

    for (uint32_t n = 0; n < bonds.size(); n++) {
      auto bond = bonds[n];

      for (uint32_t i = 0; i < system_size(); i++) {
        auto idx = tensor_idx(i);
        if (!bond.filter(idx[0], idx[1], idx[2], idx[3])) {
          continue;
        }

        std::optional<uint32_t> j_opt = get_neighbor(i, bond.d1, bond.d2, bond.d3, bond.ds);
        if (j_opt) {
          uint32_t j = j_opt.value();
          neighbors[i].push_back(Bond{j, n});
        }
      }
    }
  }

  Eigen::Vector3d position(int i1, int i2, int i3, int s) const {
    Eigen::Vector3d v = mod(i1, dx.N)*dx.v + mod(i2, dy.N)*dy.v + mod(i3, dz.N)*dz.v;
    if (s > 0) {
      v += sublattice_vectors[mod(s - 1, num_sublattices)];
    }

    return v;
  }

  Eigen::Vector3d position(int i) const {
    Eigen::Vector4i idx = tensor_idx(i);
    return position(idx[0], idx[1], idx[2], idx[3]);
  }

  uint32_t flat_idx(int n1, int n2, int n3, int s) const {
    return n1 + dx.N*(n2 + dy.N*(n3 + dz.N*s));
  }

  Eigen::Vector4i tensor_idx(int i) const {
    int n1 = i % dx.N;
    i = i / dx.N;
    int n2 = i % dy.N;
    i = i / dy.N;
    int n3 = i % dz.N;
    i = i / dz.N;
    int s = i % num_sublattices;
    Eigen::Vector4i v; v << n1, n2, n3, s;
    return v;
  }

  std::optional<uint32_t> get_neighbor(int i, int d1, int d2, int d3, int ds) const {
    Eigen::Vector4i idx = tensor_idx(i);

    int nx = idx[0] + d1;
    int ny = idx[1] + d2;
    int nz = idx[2] + d3;
    int ns = idx[3] + ds;

    if (dx.bc == BoundaryCondition::Open) {
      if (nx < 0 || nx >= dx.N) {
        return std::nullopt;
      }
    } else if (dx.bc == BoundaryCondition::Periodic) {
      nx = mod(nx, dx.N);
    }

    if (dy.bc == BoundaryCondition::Open) {
      if (ny < 0 || ny >= dy.N) {
        return std::nullopt;
      }
    } else if (dy.bc == BoundaryCondition::Periodic) {
      ny = mod(ny, dy.N);
    }

    if (dz.bc == BoundaryCondition::Open) {
      if (nz < 0 || nz >= dz.N) { 
        return std::nullopt;
      }
    } else if (dz.bc == BoundaryCondition::Periodic) {
      nz = mod(nz, dz.N);
    }

    ns = mod(ns, num_sublattices);

    return flat_idx(nx, ny, nz, ns);
  }

  double intensity(const Eigen::Vector3d& Q) const {
    constexpr int dim = Spin::RowsAtCompileTime;
    Eigen::VectorXcd structure_factor = Eigen::VectorXcd::Zero(dim);
    for (uint32_t i = 0; i < system_size(); i++) {
      structure_factor += spins[i]*std::exp(std::complex<double>(0., Q.dot(position(i))));
    }

    return std::pow(std::abs(structure_factor.norm())/system_size(), 2);
  }

  std::vector<double> intensity(const Eigen::Vector3d& q1, const Eigen::Vector3d& q2, size_t resolution) const {
    std::vector<double> intensity_samples(resolution);
    for (uint32_t i = 0; i < resolution; i++) {
      double L = i/resolution;
      Eigen::Vector3d q = L*q1 + (1 - L)*q2;
      intensity_samples[i] = intensity(q);
    }

    return intensity_samples;
  }

  std::array<std::vector<Eigen::MatrixXd>, 3> get_twist_matrices(double alpha) const {
    constexpr int dim = Spin::RowsAtCompileTime;

    std::vector<Eigen::MatrixXd> R1s;
    std::vector<Eigen::MatrixXd> R2s;
    std::vector<Eigen::MatrixXd> R3s;

    for (size_t i = 0; i < bonds.size(); i++) {
      auto bond = bonds[i];
      // TODO check this; unsigned int is likely causing issues
      Eigen::Vector3d v = position(bond.d1, bond.d2, bond.d3, bond.ds);
      double f = v[0]/v.norm() * alpha;
      Eigen::MatrixXd R = Eigen::MatrixXd::Identity(dim, dim);
      R(0, 0) =  std::cos(f);
      R(0, 1) = -std::sin(f);
      R(1, 0) =  std::sin(f);
      R(1, 1) =  std::cos(f);

      R1s.push_back(R);
      R2s.push_back(R*R);
      R3s.push_back(R*R*R);
    }

    return {R1s, R2s, R3s};
  }
};


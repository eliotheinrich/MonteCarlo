#pragma once

#include "Random.hpp"

#include <iostream>
#include <vector>
#include <iterator>
#include <set>
#include <map>
#include <algorithm>
#include <cstdint>
#include <random>
#include <climits>
#include <optional>
#include <ranges>

struct DirectedTag {};
struct UndirectedTag {};

template <typename T, bool IsVoidWeight>
struct NeighborView {};

template <typename T>
struct NeighborView<T, true> {
  const std::set<uint32_t>& s;
  auto begin() const { return s.begin(); }
  auto end()   const { return s.end(); }
};

template <typename T>
struct NeighborView<T, false> {
  const std::map<uint32_t, T>& m;
  auto begin() const { return std::views::keys(m).begin(); }
  auto end()   const { return std::views::keys(m).end(); }
};

template <typename V = int, typename T = void, typename DirTag = UndirectedTag>
class Graph {
  using EdgeContainer = std::conditional_t<
    std::is_void_v<T>,
    std::set<uint32_t>,
    std::map<uint32_t, T>
  >;

  public:
    std::vector<EdgeContainer> edges;
    std::vector<V> vals;

    uint32_t num_vertices;

    Graph() : num_vertices(0) {}

    Graph(uint32_t num_vertices) : num_vertices(0) { 
      for (uint32_t i = 0; i < num_vertices; i++) {
        add_vertex();
      }
    }

    Graph(const Graph &g) : num_vertices(0) {
      for (uint32_t i = 0; i < g.num_vertices; i++) {
        add_vertex(static_cast<V>(g.vals[i]));
      }

      for (uint32_t i = 0; i < g.num_vertices; i++) {
        if constexpr (std::is_void_v<T>) {
          for (auto const &j : g.edges[i]) {
            add_edge(i, j);
          }
        } else {
          for (auto const &[j, w] : g.edges[i]) {
            add_edge(i, j, w);
          }
        }
      }
    }

    static Graph<V, T, UndirectedTag> erdos_renyi_graph(uint32_t num_vertices, double p) {
      Graph<V, T, UndirectedTag> g(num_vertices);

      for (uint32_t i = 0; i < num_vertices; i++) {
        for (uint32_t j = i+1; j < num_vertices; j++) {
          if (randf() < p) {
            g.toggle_edge(i, j);
          }
        }
      }

      return g;
    }

    static Graph<V, T, UndirectedTag> random_regular_graph(uint32_t num_vertices, size_t k, uint32_t max_depth=0) {
      if (num_vertices*k % 2 == 1) {
        throw std::invalid_argument("To generate random regular graph, num_vertices*k must be even.");
      }

      if (k >= num_vertices) {
        throw std::invalid_argument("k must be less than num_vertices.");
      }

      Graph<V, T, UndirectedTag> buckets(num_vertices*k);
      Graph<V, T, UndirectedTag> g(num_vertices);
      std::vector<size_t> sites(num_vertices*k);

      recursive_random_regular_graph(buckets, g, sites, max_depth, 0);
      return g;
    }

    static Graph<V, T, UndirectedTag> scale_free_graph(uint32_t num_vertices, double alpha) {
      Graph<V, T, UndirectedTag> g(num_vertices);

      std::minstd_rand rng(randi());
      std::uniform_real_distribution<> dis(0.0, 1.0);

      std::vector<uint32_t> degrees(num_vertices);
      for (uint32_t i = 0; i < num_vertices; i++) {
        double u = dis(rng);
        double x = std::pow(1.0 - u * (1.0 - std::pow(1.0/num_vertices, 1.0 - alpha)), 1.0 / (1.0 - alpha));

        degrees[i] = x * (num_vertices - 1) + 1;
      }

      // Sort in reverse order
      std::sort(degrees.begin(), degrees.end(), [](uint32_t a, uint32_t b) { return a > b; });

      std::vector<uint32_t> all_vertices(num_vertices);
      std::iota(all_vertices.begin(), all_vertices.end(), 0);


      for (uint32_t i = 0; i < num_vertices-1; i++) {
        std::vector<uint32_t> random_vertices;
        uint32_t j = i+1;
        uint32_t residual_vertices = degrees[i] - g.degree(i);
        while (random_vertices.size() < residual_vertices && j < num_vertices) {
          if (g.degree(j) < degrees[j]) {
            random_vertices.push_back(j);
          }

          j++;
        }

        std::shuffle(random_vertices.begin(), random_vertices.end(), rng);

        for (uint32_t j = 0; j < random_vertices.size(); j++) {
          g.add_edge(i, random_vertices[j]);
        }
      }

      return g;
    }

    std::string to_string() const {
      std::string s = "";
      for (uint32_t i = 0; i < num_vertices; i++) {
        s += std::format("[{}] {} -> ", vals[i], i); 
        if constexpr (std::is_void_v<T>) {
          for (auto const& v : edges[i]) {
            s += std::format("({}) ", v);
          }
        } else {
          for (auto const&[v, w] : edges[i]) {
            s += std::format("({}: {}) ", v, w);
          }
        }
        if (i != num_vertices - 1) {
          s += "\n";
        }
      }
      return s;
    }

    size_t size() const {
      return num_vertices;
    }

    void add_vertex(std::optional<V> val_opt=std::nullopt) {
      num_vertices++;
      if constexpr (std::is_void_v<T>) {
        edges.push_back(std::set<uint32_t>());
      } else {
        edges.push_back(std::map<uint32_t, T>());
      }
      vals.push_back(val_opt.value_or(V()));
    }

    void remove_vertex(uint32_t u) {
      if (u >= num_vertices) {
        throw std::runtime_error(std::format("Cannot remove vertex {} from a graph with {} vertices.", u, num_vertices));
      }
      num_vertices--;
      edges.erase(edges.begin() + u);
      vals.erase(vals.begin() + u);
      for (uint32_t i = 0; i < num_vertices; i++) {
        edges[i].erase(u);
      }

      for (uint32_t i = 0; i < num_vertices; i++) {
        if constexpr (std::is_void_v<T>) {
          std::set<uint32_t> new_edges;

          for (auto const &j : edges[i]) {
            if (j > u) {
              new_edges.insert(j - 1);
            } else {
              new_edges.insert(j);
            }
          }
          edges[i] = std::move(new_edges);
        } else {
          std::map<uint32_t, T> new_edges;
          for (auto const &[j, w] : edges[i]) {
            if (j > u) {
              new_edges.emplace(j-1, w);
            } else {
              new_edges.emplace(j, w);
            }
          }
          edges[i] = std::move(new_edges);
        }
      }
    }

    void set_val(uint32_t u, const V& val) {
      if (u >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertex {} for graph with {} vertices.", u, num_vertices));
      }

      vals[u] = val;
    }

    const V& get_val(uint32_t u) const { 
      if (u >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertex {} for graph with {} vertices.", u, num_vertices));
      }
      return vals[u]; 
    }

    auto edges_of(uint32_t u) const {
      if (u >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertex {} for graph with {} vertices.", u, num_vertices));
      }

      if constexpr (std::is_void_v<T>) {
        return NeighborView<T, true>{edges[u]};
      } else {
        return NeighborView<T, false>{edges[u]};
      }
    }

    std::vector<uint32_t> neighbors(uint32_t u) const {
      std::vector<uint32_t> neighbors(edges_of(u).begin(), edges_of(u).end());
      std::sort(neighbors.begin(), neighbors.end());
      return neighbors;
    }

    bool contains_edge(uint32_t u1, uint32_t u2) const {
      if (u1 >= num_vertices || u2 >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertices {} and {} for graph with {} vertices.", u1, u2, num_vertices));
      }
      if constexpr (std::is_same_v<DirTag, DirectedTag>) {
        return edges[u1].contains(u2);
      } else {
        return edges[u1].contains(u2) && edges[u2].contains(u1);
      }
    }

    T edge_weight(uint32_t u1, uint32_t u2) const {
      static_assert(!std::is_void_v<T>, "Unweighted graph has no edge_weights.");
      if (u1 >= num_vertices || u2 >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertices {} and {} for graph with {} vertices.", u1, u2, num_vertices));
      }
      return edges[u1].at(u2);
    }

    void set_edge_weight(uint32_t u1, uint32_t u2, std::conditional_t<std::is_void_v<T>, std::nullptr_t, T> w={}) {
      static_assert(!std::is_void_v<T>, "Unweighted graph has no edge_weights.");
      if (u1 >= num_vertices || u2 >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertices {} and {} for graph with {} vertices.", u1, u2, num_vertices));
      }
      edges[u1][u2] = w;
    }

    void add_edge(uint32_t u1, uint32_t u2, std::conditional_t<std::is_void_v<T>, std::nullptr_t, T> weight={}) {
      if (u1 >= num_vertices || u2 >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertices {} and {} for graph with {} vertices.", u1, u2, num_vertices));
      }

      if (contains_edge(u1, u2)) {
        return; 
      }

      if constexpr (std::is_void_v<T>) {
        edges[u1].insert(u2);
      } else {
        edges[u1].emplace(u2, weight);
      }

      if constexpr(std::is_same_v<DirTag, UndirectedTag>) {
        if (u1 != u2) {
          add_edge(u2, u1, weight);
        }
      }
    }

    void remove_edge(uint32_t u1, uint32_t u2) {
      if (u1 >= num_vertices || u2 >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertices {} and {} for graph with {} vertices.", u1, u2, num_vertices));
      }

      edges[u1].erase(u2);
      if constexpr (std::is_same_v<DirTag, UndirectedTag>) {
        if (u1 != u2) {
          edges[u2].erase(u1);
        }
      }
    }

    void toggle_edge(uint32_t u1, uint32_t u2) {
      if (u1 >= num_vertices || u2 >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertices {} and {} for graph with {} vertices.", u1, u2, num_vertices));
      }

      if (contains_edge(u1, u2)) {
        remove_edge(u1, u2);
      } else {
        add_edge(u1, u2);
      }
    }

    std::conditional_t<std::is_void_v<T>, bool, std::optional<T>> adjacency_matrix(uint32_t u1, uint32_t u2) const {
      if (u1 >= num_vertices || u2 >= num_vertices) {
        throw std::runtime_error(std::format("Invalid vertices {} and {} for graph with {} vertices.", u1, u2, num_vertices));
      }

      if (contains_edge(u1, u2)) {
        if constexpr (std::is_void_v<T>) {
          return true;
        } else {
          return edges[u1].at(u2);
        }
      } else {
        if constexpr (std::is_void_v<T>) {
          return false;
        } else {
          return std::nullopt;
        }
      }
    }

    uint32_t degree(uint32_t u) const {
      return edges[u].size();
    }

    void local_complement(uint32_t u) {
      for (auto const& v1 : edges_of(u)) {
        for (auto const& v2 : edges_of(u)) {
          if (v1 < v2) {
            toggle_edge(v1, v2);
          }
        }
      }
    }

    size_t num_edges() const {
      size_t n = 0;
      for (size_t i = 0; i < num_vertices; i++) {
        n += edges[i].size();
      }
      return n;
    }

    Graph<V, T, DirTag> subgraph(const std::vector<uint32_t>& sites) const {
      Graph<V, T, DirTag> g(sites.size());
      for (size_t i = 0; i < sites.size(); i++) {
        size_t a = sites[i];
        if constexpr (std::is_void_v<T>) {
          for (auto const& b: edges[a]) {
            for (size_t j = 0; j < sites.size(); j++) {
              if (sites[j] == b) {
                g.add_edge(i, j);
                break;
              }
            }
          }
        } else {
          for (auto const& [b, w] : edges[a]) {
            for (size_t j = 0; j < sites.size(); j++) {
              if (sites[j] == b) {
                g.add_edge(i, j, w);
                break;
              }
            }
          }
        }
      }

      return g;
    }

    Graph<bool, T, DirectedTag> partition(const std::vector<uint32_t> &nodes) const {
      std::set<uint32_t> nodess;
      std::copy(nodes.begin(), nodes.end(), std::inserter(nodess, nodess.end()));
      Graph<bool, T, DirectedTag> new_graph;
      std::map<uint32_t, uint32_t> new_vertices;

      for (const uint32_t a : nodess) {
        if (!degree(a)) {
          continue;
        }

        new_vertices.emplace(a, new_vertices.size());
        new_graph.add_vertex(true);
        for (auto const &b : edges_of(a)) {
          if (nodess.count(b)) {
            continue;
          }

          if (!new_vertices.count(b)) {
            new_vertices.emplace(b, new_vertices.size());
            new_graph.add_vertex(false);
          }

          new_graph.add_edge(new_vertices[a], new_vertices[b]);
        }
      }

      bool continue_deleting = true;

      while (continue_deleting) {
        continue_deleting = false;
        for (uint32_t i = 0; i < new_graph.num_vertices; i++) {
          if (!new_graph.degree(i)) {
            new_graph.remove_vertex(i);
            continue_deleting = true;
            break;
          }
        }
      }

      return new_graph;
    }

    size_t num_loops() const {
      auto components = component_partition();

      size_t n = 0;
      for (auto const &component : components) {
        std::vector<uint32_t> sites(component.begin(), component.end());
        auto graph = subgraph(sites);
        n += graph.num_edges()/2 - graph.num_vertices + 1;
      }
      
      return n;
    }

    std::pair<bool, std::vector<uint32_t>> path(uint32_t s, uint32_t t) const requires (std::is_integral_v<T>) {
      std::vector<uint32_t> stack;
      stack.push_back(s);

      std::set<uint32_t> visited;
      std::map<uint32_t, uint32_t> parent;

      while (!stack.empty()) {
        uint32_t v = *(stack.end()-1);
        stack.pop_back();
        if (visited.count(v)) {
          continue;
        }

        visited.insert(v);
        for (auto const &w : edges_of(v)) {
          if (edge_weight(v, w) > 0) {
            if (!visited.count(w)) {
              parent.emplace(w, v);
            }
            if (w == t) {
              // Done; re-use stack
              stack.clear();
              uint32_t node = t;
              while (parent.count(node)) {
                stack.push_back(node);
                node = parent[node];
              }
              stack.push_back(s);
              std::reverse(stack.begin(), stack.end());

              return std::pair(true, stack);
            }
            stack.push_back(w);
          }
        }
      }
      return std::pair(false, std::vector<uint32_t>());
    }

    T max_flow(std::vector<uint32_t> &sources, std::vector<uint32_t> &sinks) const requires (std::is_integral_v<T>) {
      Graph g(*this);

      g.add_vertex();
      uint32_t s = g.num_vertices - 1;
      g.add_vertex();
      uint32_t t = g.num_vertices - 1;

      for (const auto& i : sources) {
        g.add_edge(s, i, INT_MAX);
      }
      for (const auto& i : sinks) {
        g.add_edge(i, t, INT_MAX);
      }

      T flow = g.max_flow(s, t);

      return flow;
    }

    T max_flow(uint32_t s, uint32_t t) const requires (std::is_integral_v<T>) {
      Graph residual_graph(*this);

      for (uint32_t i = 0; i < num_vertices; i++) {
        for (auto const &w : residual_graph.edges_of(i)) {
          residual_graph.add_edge(w, i, 0);
        }
      }

      std::pair<bool, std::vector<uint32_t>> p = residual_graph.path(s, t);
      bool path_exists = p.first;
      std::vector<uint32_t> path_nodes = p.second;

      while (path_exists) {
        T min_weight = INT_MAX;
        for (uint32_t j = 0; j < path_nodes.size() - 1; j++) {
          T weight = residual_graph.edge_weight(path_nodes[j], path_nodes[j+1]);
          if (weight < min_weight) {
            min_weight = weight;
          }
        }

        for (uint32_t j = 0; j < path_nodes.size() - 1; j++) {
          uint32_t u = path_nodes[j];
          uint32_t v = path_nodes[j+1];
          residual_graph.set_edge_weight(v, u, residual_graph.edge_weight(v, u) + min_weight);
          residual_graph.set_edge_weight(u, v, residual_graph.edge_weight(u, v) - min_weight);
        }

        p = residual_graph.path(s, t);
        path_exists = p.first;
        path_nodes = p.second;
      }

      T flow = 0;
      for (auto const &[v, w] : residual_graph.edges[s]) {
        flow += edge_weight(s, v) - w;
      }

      return flow;
    }

    std::set<uint32_t> component(uint32_t i) const {
      std::vector<uint32_t> stack;
      stack.push_back(i);
      std::set<uint32_t> visited;

      while (!stack.empty()) {
        uint32_t v = *(stack.end()-1);
        stack.pop_back();
        if (!visited.count(v)) {
          visited.insert(v);
          for (auto const &w : edges_of(v)) {
            stack.push_back(w);
          }
        }
      }

      return visited;
    }


    // Graph properties
    std::vector<uint32_t> compute_degree_counts() const {
      std::vector<uint32_t> counts(num_vertices, 0);
      for (uint32_t i = 0; i < num_vertices; i++) {
        counts[degree(i)]++;
      }
      return counts;
    }

    std::vector<uint32_t> compute_neighbor_degree_counts() const {
      std::vector<uint32_t> counts(num_vertices, 0);
      for (uint32_t i = 0; i < num_vertices; i++) {
        if (degree(i) > 0) {
          uint32_t v = randi() % degree(i);
          counts[degree(v)]++;
        }
      }
      return counts;
    }

    double average_component_size() const {
      std::set<uint32_t> to_check;
      for (uint32_t i = 0; i < num_vertices; i++) {
        to_check.insert(i);
      }
      double avg = 0.;

      while (!to_check.empty()) {
        uint32_t i = *to_check.begin(); // pop
        to_check.erase(to_check.begin());

        auto connected_component = component(i);
        uint32_t component_size = connected_component.size();
        for (auto v : connected_component) {
          if (to_check.count(v)) {
            to_check.erase(v);
          }
        }

        avg += component_size*component_size;
      }

      return avg/num_vertices;
    }

    uint32_t max_component_size() const {
      std::set<uint32_t> to_check;
      for (uint32_t i = 0; i < num_vertices; i++) {
        to_check.insert(i);
      }

      uint32_t max_cluster_size = 0;

      while (!to_check.empty()) {
        uint32_t i = *to_check.begin(); // pop
        to_check.erase(to_check.begin());

        auto connected_component = component(i);
        uint32_t component_size = connected_component.size();
        if (component_size > max_cluster_size) {
          max_cluster_size = component_size;
        }

        for (auto v : connected_component) {
          if (to_check.count(v)) {
            to_check.erase(v);
          }
        }
      }

      return max_cluster_size;
    }

    double local_clustering_coefficient(uint32_t i) const requires (std::is_arithmetic_v<T> || std::is_void_v<T>) {
      uint32_t ki = degree(i);
      if (ki == 0 || ki == 1) {
        return 0.;
      }

      uint32_t c = 0;
      for (uint32_t j = 0; j < num_vertices; j++) {
        for (uint32_t k = 0; k < num_vertices; k++) {
          if constexpr (std::is_void_v<T>) {
            c += adjacency_matrix(i, j)*adjacency_matrix(j, k)*adjacency_matrix(k, i);
          } else {
            std::optional<T> c1 = adjacency_matrix(i, j);
            std::optional<T> c2 = adjacency_matrix(j, k);
            std::optional<T> c3 = adjacency_matrix(k, i);
            if (c1 && c2 && c3) {
              c += c1.value() * c2.value() * c3.value();
            }
          }
        }
      }

      return c/(ki*(ki - 1));
    }
      
    double global_clustering_coefficient() const requires (std::is_arithmetic_v<T> || std::is_void_v<T>) {
      double c = 0.;
      for (uint32_t i = 0; i < num_vertices; i++) {
        c += local_clustering_coefficient(i);
      }
      return c/num_vertices;
    }

    std::vector<std::set<uint32_t>> component_partition() const {
      std::vector<std::set<uint32_t>> components;

      std::set<uint32_t> to_check;
      for (uint32_t i = 0; i < num_vertices; i++) {
        to_check.insert(i);
      }

      while (!to_check.empty()) {
        uint32_t i = *to_check.begin(); // pop
        to_check.erase(to_check.begin());

        auto new_component = component(i);
        components.push_back(new_component);

        for (auto v : new_component) {
          if (to_check.count(v)) {
            to_check.erase(v);
          }
        }
      }

      return components;
    }

    double percolation_probability() const {
      std::set<uint32_t> to_check;
      for (uint32_t i = 0; i < num_vertices; i++) {
        to_check.insert(i);
      }

      uint32_t max_cluster_size = 0;

      while (!to_check.empty()) {
        uint32_t i = *to_check.begin(); // pop
        to_check.erase(to_check.begin());

        auto connected_component = component(i);
        uint32_t component_size = connected_component.size();
        if (component_size > max_cluster_size) {
          max_cluster_size = component_size;
        }

        for (const auto& v : connected_component) {
          if (to_check.count(v)) {
            to_check.erase(v);
          }
        }
      }

      return double(max_cluster_size)/num_vertices;
    }

    private:
      static void recursive_random_regular_graph(
        Graph<V, T, DirTag>& buckets, 
        Graph<V, T, DirTag>& g, 
        std::vector<size_t>& sites, 
        uint32_t max_depth, 
        uint32_t depth
      ) {
        if (max_depth != 0) {
          if (depth > max_depth) {
            throw std::invalid_argument("Maximum depth reached.");
          }
        }

        // Create pairs
        buckets = Graph<V, T, DirTag>(buckets.num_vertices);
        g = Graph<V, T, DirTag>(g.num_vertices);
        std::iota(sites.begin(), sites.end(), 0);
        std::minstd_rand rng(randi());
        std::shuffle(sites.begin(), sites.end(), rng);

        for (size_t i = 0; i < sites.size()/2; i++) {
          buckets.add_edge(sites[i], sites[i + sites.size()/2]);
        }

        size_t num_vertices = g.num_vertices;
        size_t k = buckets.num_vertices/num_vertices;

        // Collapse nodes
        for (size_t v1 = 0; v1 < num_vertices*k; v1++) {
          for (auto const& [v2, _] : buckets.edges[v1]) {
            // Only consider each edge once
            if (v2 < v1) {
              continue;
            }

            size_t i = v1 / k;
            size_t j = v2 / k;

            if (i == j || g.contains_edge(i, j)) {
              recursive_random_regular_graph(buckets, g, sites, max_depth, depth + 1);
              return;
            }

            g.add_edge(i, j);
          }
        }
      }
};

template <typename V, typename T=void>
using DirectedGraph = Graph<V, T, DirectedTag>;

template <typename V, typename T=void>
using UndirectedGraph = Graph<V, T, UndirectedTag>;


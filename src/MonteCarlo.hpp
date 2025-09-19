#pragma once

#include <cmath>
#include <map>
#include <string>
#include <variant>
#include <iostream>

#include "Random.hpp"

using Parameter = std::variant<std::string, int, double>;
using Params = std::map<std::string, Parameter>;

template <class T>
T get(Params &params, const std::string& key, T defaultv) {
  if (params.count(key)) {
    return std::get<T>(params[key]);
  }

  params[key] = Parameter{defaultv};
  return defaultv;
}

template <class T>
T get(const Params &params, const std::string& key) {
  if (!params.count(key)) {
    throw std::runtime_error(std::format("Key \"{}\" not found in Params.", key));
  }
  return std::get<T>(params.at(key));
}

static inline uint32_t mod(int a, int b) {
  int c = a % b;
  return (c < 0) ? c + b : c;
}

enum BoundaryCondition { Periodic, Open };

typedef std::pair<uint32_t, int> Bond;

inline BoundaryCondition parse_boundary_condition(std::string s) {
  if (s == "periodic") {
    return BoundaryCondition::Periodic;
  } else if (s == "open") {
    return BoundaryCondition::Open;
  } else {
    std::string error_message = "Invalid boundary condition: " + s + "\n";
    throw std::invalid_argument(error_message);
  }
}

struct Color {
  float r;
  float g;
  float b;
  float w;
};

template <>
struct std::formatter<Color> : std::formatter<std::string> {
  auto format(const Color& c, std::format_context& ctx) const {
    return std::formatter<std::string>::format(
      std::format("{}, {}, {}, {}", c.r, c.g, c.b, c.w), ctx);
  }
};

class Texture {
  public:
    size_t n;
    size_t m;
    std::vector<float> texture;

    Texture(const std::vector<float>& data, size_t n, size_t m) : n(n), m(m), texture(data) {
      if (len() != data.size()) {
        throw std::runtime_error("Invalid texture dimensions passed to Texture.");
      }
    }

    Texture(size_t n, size_t m) : n(n), m(m), texture(4*n*m) {}

    Texture(const std::vector<std::vector<std::vector<float>>>& data) {
      n = data.size();
      if (n > 0) {
        m = data[0].size();
        for (size_t k = 0; k < n; k++) {
          if (data[k].size() != m) {
            throw std::runtime_error("Data not square.");
          }
          for (size_t l = 0; l < m; l++) {
            if (data[k][l].size() != 4) {
              throw std::runtime_error("Data is not RGBA formatted.");
            }
          }
        }
      }

      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
          set(i, j, {data[i][j][0], data[i][j][1], data[i][j][2], data[i][j][3]});
        }
      }
    }

    Texture()=default;

    void set(size_t i, size_t j, Color color) {
      size_t idx = 4*(i + j*m);
      texture[idx]   = color.r;
      texture[idx+1] = color.g;
      texture[idx+2] = color.b;
      texture[idx+3] = color.w;
    }
    
    size_t len() const {
      return 4*n*m;
    }

    const float* data() const {
      return texture.data();
    }
};

class MonteCarloSimulator {
  // Most basic Monte-Carlo model to be simulated must have some notion of energy
  // as well as a mutation data structure. Specifics must be supplied by child classes.
  public:
    MonteCarloSimulator()=default;

    MonteCarloSimulator(Params &params, uint32_t num_threads) : num_threads(num_threads) {
      temperature = get<double>(params, "temperature");

      if (temperature < 0) {
        throw std::runtime_error(std::format("Provided temperature {:.3f} is negative.", temperature));
      }
    }

    virtual ~MonteCarloSimulator();

    // Implement Simulator methods but introduce MonteCarlo methods
    virtual void timesteps(uint32_t num_steps) {
      uint64_t num_updates = system_size()*num_steps;
      for (uint64_t i = 0; i < num_updates; i++) {
        generate_mutation();
        double dE = energy_change();

        double rf = randf();
        if (rf < std::exp(-dE/temperature)) {
          accept_mutation();
        } else {
          reject_mutation();
        }
      }
    }

    virtual void key_callback(int key);

    // To be overridden by child classes
    virtual double energy() const=0;
    virtual double energy_change()=0;
    virtual void generate_mutation()=0;
    virtual void accept_mutation()=0;
    virtual void reject_mutation()=0;
    virtual uint64_t system_size() const=0;

    virtual Texture get_texture() const {
      throw std::runtime_error("Called get_texture on a C++ simulator that does not implement it.");
    }

  protected:
    double temperature;

  private:
    uint32_t num_threads;
};


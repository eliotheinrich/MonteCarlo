#pragma once

#include <random>
#include <stdexcept>
#include <format>

class Random {
  public:
    static Random& get_instance() {
      thread_local static Random instance;
      return instance;
    }
  private:
    uint32_t seed;
    std::minstd_rand rng;
    Random() {
      thread_local std::random_device gen;
      seed = gen();
      rng.seed(seed);
    } 

  public:
    Random(const Random&) = delete;
    Random& operator=(const Random&) = delete;

    static void seed_rng(uint32_t s) {
      Random& instance = get_instance();
      instance.seed = s;
      instance.rng.seed(s);
    }

    static uint32_t get_seed() {
      return Random::get_instance().seed;
    }

    uint32_t rand() {
      return rng();
    }
};

inline static uint32_t randi() {
  return Random::get_instance().rand();
}

inline static uint32_t randi(uint32_t min, uint32_t max) {
  if (max <= min) {
    throw std::runtime_error(std::format("Max must be greater then min for generating random integers."));
  }
  return randi() % (max - min) + min;
}

inline static double randf() {
  return static_cast<double>(randi())/static_cast<double>(RAND_MAX);
}

inline static double randf(double min, double max) {
  if (max <= min) {
    throw std::runtime_error(std::format("Max must be greater then min for generating random floats."));
  }
  return randf() * (max - min) + min;
}

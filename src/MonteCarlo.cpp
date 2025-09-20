#include "MonteCarlo.hpp"

#include <GLFW/glfw3.h>

void MonteCarloSimulator::key_callback(const std::string& key) {
  if (key == "ArrowUp") {
    std::cout << std::format("Lowering temperature: T = {:.2f}\n", temperature);
    temperature = std::max(0.0, temperature - 0.1);
  } else if (key == "ArrowDown") {
    std::cout << std::format("Raising temperature: T = {:.2f}\n", temperature);
    temperature = temperature + 0.1;
  }
}

MonteCarloSimulator::~MonteCarloSimulator()=default;

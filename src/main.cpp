#include "Ising/SquareIsingModel.h"
#include "MonteCarlo.hpp"
#include "Spin2d/SquareXYModel.h"
#include "SimulatorDisplay/Animator.h"
#include <chrono>


std::shared_ptr<MonteCarloSimulator> prepare_model() {
  Params params;
  params.emplace("system_size", 256);
  params.emplace("J", 1.0);
  params.emplace("B", 0.0);
  params.emplace("Bp", 0.0);
  params.emplace("temperature", 1.0);
  params.emplace("cluster_update", 0);

  return std::make_shared<SquareXYModel>(params, 1);
}

int main(int argc, char** argv) {
  Animator animator(prepare_model());
  animator.start(900, 900);

  return 0;
}

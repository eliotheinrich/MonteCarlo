#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Frame.h>

#include "Spin3d/Selenium.h"
#include "Spin2d/SquareXYModel.h"
#include "Spin3d/XXZHeis.h"

using namespace dataframe;
using namespace dataframe::utils;

template <class SimulatorType>
void emulate_simulation(ExperimentParams& params, size_t num_threads=1) {
  SimulatorType sim(params, num_threads);

  size_t equilibration_timesteps = get<int>(params, "equilibration_timesteps", 0);
  size_t sampling_timesteps = get<int>(params, "sampling_timesteps", 0);
  size_t measurement_freq = get<int>(params, "measurement_freq", 1);
  size_t temporal_avg = get<int>(params, "temporal_avg", true);

  size_t num_steps = sampling_timesteps / measurement_freq;

  sim.timesteps(equilibration_timesteps);

  for (size_t i = 0; i < num_steps; i++) {
    sim.timesteps(measurement_freq);
    auto samples = sim.take_samples();
  }

  return;
}

bool test_selenium() {
  ExperimentParams params;
  params["system_size"] = static_cast<int>(3);
  params["temperature"] = 1.0;

  params["K"]   = -0.9;
  params["J1"]  =  1.0;
  params["J2"]  = -0.6;
  params["J2p"] = -0.6;
  params["D1"]  = -0.1;
  params["D2"]  = -0.1;

  params["equilibration_timesteps"] = static_cast<int>(100);
  params["cluster_update"]          = static_cast<int>(0);
  params["sampling_timesteps"]      = static_cast<int>(1000);
  params["temporal_avg"]            = static_cast<int>(0);
  params["measurement_freq"]        = static_cast<int>(2);

  params["sample_magnetization"]            = static_cast<int>(0);
  params["sample_energy"]                   = static_cast<int>(0);
  params["sample_sublattice_magnetization"] = static_cast<int>(0);

  params["sample_structure_factor"] = static_cast<int>(1);


  emulate_simulation<SeleniumModel>(params, 1);

  return true;
}

bool test_xxz() {
  ExperimentParams params;
  params["system_size"] = static_cast<int>(4);
  params["temperature"] = 0.2;
  params["initial_temperature"] = 2.0;

  params["J"] = -1.0;
  params["A"] = -0.8;
  params["K"] = 0.0;

  params["equilibration_timesteps"] = static_cast<int>(1000);
  params["cluster_update"]          = static_cast<int>(0);
  params["sampling_timesteps"]      = static_cast<int>(10000);
  params["temporal_avg"]            = static_cast<int>(1);
  params["measurement_freq"]        = static_cast<int>(5);

  params["sample_magnetization"]    = static_cast<int>(0);
  params["sample_energy"]           = static_cast<int>(0);
  params["sample_helicity"]         = static_cast<int>(0);
  params["sample_structure_factor"] = static_cast<int>(1);

  emulate_simulation<XXZHeis>(params, 1);

  return true;
}


int main() {
  test_selenium();
  //test_xxz();
}

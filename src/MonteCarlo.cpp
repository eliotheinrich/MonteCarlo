#include "MonteCarlo.h"
#include <iostream>

MonteCarloSimulator::MonteCarloSimulator(Params &params, MCModel *model) {
    temperature = params.getf("temperature");
    init_temperature = params.getf("initial_temperature", temperature);
    num_cooling_steps = params.getf("num_cooling_steps", DEFAULT_NUM_COOLING_STEPS);
    cooling_schedule = parse_cooling_schedule(params.gets("cooling_schedule", DEFAULT_COOLING_SCHEDULE));
    random_seed = params.geti("random_seed", DEFAULT_RANDOM_SEED);

    if (random_seed == -1) rng = new std::minstd_rand();
    else rng = new std::minstd_rand(random_seed);

    this->model = model->clone(params);
}

void MonteCarloSimulator::init_state() {
    model->init();
}

void MonteCarloSimulator::timesteps(uint num_steps) {
    for (uint i = 0; i < num_steps; i++) {
        model->generate_mutation();
        double dE = model->energy_change();

        double rf = randf();
        if (rf < std::exp(-dE/model->temperature)) {
            model->accept_mutation();
        } else {
            model->reject_mutation();
        }
    }
}

void MonteCarloSimulator::equilibration_timesteps(uint num_steps) {
    for (uint i = 0; i < num_cooling_steps; i++) {
        timesteps(num_steps/num_cooling_steps);
        switch (cooling_schedule) {
            case(CoolingSchedule::Constant) : model->temperature = const_T(i, num_cooling_steps, init_temperature, temperature);
            case(CoolingSchedule::Linear) : model->temperature = linear_T(i, num_cooling_steps, init_temperature, temperature);
            case(CoolingSchedule::Trig) : model->temperature = trig_T(i, num_cooling_steps, init_temperature, temperature);
        }
    }
    model->temperature = temperature;
}
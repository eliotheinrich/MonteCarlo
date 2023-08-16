#include "MonteCarlo.h"
#include <iostream>

#define DEFAULT_COOLING_SCHEDULE "constant"
#define DEFAULT_NUM_COOLING_UPDATES 100
#define DEFAULT_RANDOM_SEED -1

MonteCarloSimulator::MonteCarloSimulator(Params &params, std::shared_ptr<MCModel> model) : Simulator(params) {
    temperature = get<double>(params, "temperature");
    init_temperature = get<double>(params, "initial_temperature", temperature);
    num_cooling_updates = get<int>(params, "num_cooling_updates", DEFAULT_NUM_COOLING_UPDATES);
    cooling_schedule = parse_cooling_schedule(get<std::string>(params, "cooling_schedule", DEFAULT_COOLING_SCHEDULE));

    this->model = std::move(model);
}

void MonteCarloSimulator::init_state() {
    model->init();
}

void MonteCarloSimulator::timesteps(uint num_steps) {
    ull num_updates = model->system_size()*num_steps;
    for (ull i = 0; i < num_updates; i++) {
        model->generate_mutation();
        double dE = model->energy_change();

        double rf = randf();
        if (rf < std::exp(-dE/model->temperature))
            model->accept_mutation();
        else
            model->reject_mutation();
    }
}

void MonteCarloSimulator::equilibration_timesteps(uint num_steps) {
    uint steps_per_update = num_steps / num_cooling_updates;
    model->temperature = init_temperature;
    for (uint i = 0; i < num_cooling_updates; i++) {
        timesteps(steps_per_update);
        switch (cooling_schedule) {
            case(CoolingSchedule::Constant) : model->temperature = const_T(i, num_cooling_updates, init_temperature, temperature);
            case(CoolingSchedule::Linear) : model->temperature = linear_T(i, num_cooling_updates, init_temperature, temperature);
            case(CoolingSchedule::Trig) : model->temperature = trig_T(i, num_cooling_updates, init_temperature, temperature);
        }
    }

    model->temperature = temperature;
    timesteps(num_steps % num_cooling_updates);
}
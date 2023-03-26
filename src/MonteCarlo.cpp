#include "MonteCarlo.h"
#include <iostream>

MonteCarloSimulator::MonteCarloSimulator(Params &params, std::unique_ptr<MCModel> model) : Simulator(params) {
    temperature = params.get<float>("temperature");
    init_temperature = params.get<float>("initial_temperature", temperature);
    num_cooling_steps = params.get<float>("num_cooling_steps", DEFAULT_NUM_COOLING_STEPS);
    cooling_schedule = parse_cooling_schedule(params.get<std::string>("cooling_schedule", DEFAULT_COOLING_SCHEDULE));

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
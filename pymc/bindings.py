from pymc.montecarlo_bindings import *

from dataframe import TimeConfig

simulators = {
    "simple_graph": SimpleGraphModel,
    "square_ising": SquareIsingModel,
    "square_xy": SquareXYModel,
    "trigonal_xy": TrigonalXYModel,
    "xxz_heis": XXZHeis,
    "trigonal": TrigonalModel,
}

config_types = {
}

def prepare_config(params):
    circuit_type = params["sim_type"]
    if circuit_type in simulators:
        simulator_generator = simulators[circuit_type]
        return TimeConfig(params, simulator_generator)
    else:
        config_generator = config_types[circuit_type]
        return config_generator(params)


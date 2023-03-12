#include "MonteCarlo.h"
#include <nlohmann/json.hpp>
#include <DataFrame.hpp>

#include "SquareIsingModel.h"
#include "SquareClockModel.h"
#include "SquareXYModel.h"
#include "TrigonalXYModel.h"
#include "TrigonalModel.h"
#include "XXZHeis.h"
#include "SimpleGraphModel.h"


using json = nlohmann::json;

void defaultf() {
    std::cout << "Default behavior\n";
}

bool file_valid(std::string filename) {
    uint strlen = filename.length();
    if (strlen < 6) { return false; }
    
    std::string extension = filename.substr(strlen - 5, strlen);
    std::string json_ext = ".json";
    if (extension != json_ext) { return false; }

    std::ifstream f(filename);

    return f.good();
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        defaultf();
        return 1;
    }

    if (argc != 3) std::cout << "Incorrect arguments.\n";

    std::string filename = argv[1];
    uint num_threads = std::stoi(argv[2]);
    bool valid = file_valid(filename);
    if (!valid) {
        std::cout << "Cannot find " << filename << "; aborting.\n";
        return 1;
    }

    std::ifstream f(filename);
    json data = json::parse(f);
    std::string model_type = data["model_type"];
    std::string data_filename = data["filename"];

    std::cout << "Starting job\n";

    auto params = Params::load_json(data, true);
    std::vector<Config*> configs;

    std::string data_prefix = "../data/";

    for (auto param : params) {
        TimeConfig *config = new TimeConfig(param);

        MCModel *model;
        if (model_type == "square_ising") model = new SquareIsingModel(param);
        else if (model_type == "square_xy") model = new SquareXYModel(param);
        else if (model_type == "trigonal_xy") model = new TrigonalXYModel(param);
        else if (model_type == "square_xxz") model = new XXZHeis(param);
        else if (model_type == "trigonal_heisenberg") model = new TrigonalModel(param);
        else if (model_type == "simple_graph") model = new SimpleGraphModel(param);
        // Clock models not currently supported

        Simulator *sim = new MonteCarloSimulator(param, model);
        config->init_simulator(sim->clone(param));
        configs.push_back(config);
    }

    ParallelCompute pc(configs);
    DataFrame df = pc.compute(num_threads);
    df.write_json(data_prefix + data_filename);
    std::cout << "Finishing job\n";
}
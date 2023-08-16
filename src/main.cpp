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
#include <memory>

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

    auto params = load_json(data, true);
    std::vector<std::shared_ptr<Config>> configs;

    std::string data_prefix = "../data/";

    for (auto param : params) {
        std::shared_ptr<TimeConfig> config(new TimeConfig(param));

        std::shared_ptr<MCModel> model;
        if (model_type == "square_ising") model = std::make_shared<SquareIsingModel>(param);
        else if (model_type == "square_xy") model = std::make_shared<SquareXYModel>(param); 
        else if (model_type == "trigonal_xy") model = std::make_shared<TrigonalXYModel>(param); 
        else if (model_type == "square_xxz") model = std::make_shared<XXZHeis>(param); 
        else if (model_type == "trigonal_heisenberg") model = std::make_shared<TrigonalModel>(param);
        else if (model_type == "simple_graph") model = std::make_shared<SimpleGraphModel>(param);
        // Clock models not currently supported

        std::shared_ptr<Simulator> sim = std::shared_ptr<Simulator>(new MonteCarloSimulator(param, std::move(model)));
        config->init_simulator(sim);
        configs.push_back(config);
    }

    ParallelCompute pc(configs, num_threads);
    pc.compute(true);
    pc.write_json(data_prefix + data_filename);
    std::cout << "Finishing job\n";
}
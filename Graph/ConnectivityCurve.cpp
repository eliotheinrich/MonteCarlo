#include "SimpleGraphModel.cpp"
#include "../Routines.cpp"
#include <iostream>

std::vector<double> measure_connectivity(SimpleGraphModel *model) {
    return std::vector<double>{ model->energy(), model->get_connectivity() };
}

int main(int argc, char* argv[]) {    
    std::cout << "Starting\n";
    std::string filename = argv[1];
    srand((unsigned)time( NULL ));

    const int N = 64;
	const float J = 1.;

    const int MCStep = N;


    SimpleGraphModel *model = new SimpleGraphModel(N, J);

    int resolution = 30;
    int num_samples = 50;
    int steps_per_sample = 50;


    std::vector<double> T(0);
    float Tmin = 0.1;
    float Tmax = 3.0;
    for (int i = 0; i < resolution; i++) {
        T.push_back((i)/resolution*Tmax + double(resolution - i)/resolution*Tmin);
    }

    auto data = sample_r(measure_connectivity, model, T, 100*MCStep, num_samples, steps_per_sample, 1);
    std::string header = std::to_string(N);
    write_data(&data, &T, "data.txt", header);
}
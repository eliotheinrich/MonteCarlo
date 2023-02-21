#include "SquareIsingModel.h"
#include "MonteCarlo.h"
#include <iostream>

std::map<std::string, double> sampler(SquareIsingModel *model) {
	std::map<std::string, double> data;
	data.emplace("energy", model->energy());
	data.emplace("magnetization", model->get_magnetization());
	return data;
}


int main() {

	int num_models = 30;
	std::vector<SquareIsingModel*> models;
	double Tmin = 0.1;
	double Tmax = 5.0;

	for (int i = 0; i < num_models; i++) {
		models.push_back(new SquareIsingModel(10, 1, 1., 0.));
		models[i]->T = Tmin*i/double(num_models) + Tmax*(num_models - i)/double(num_models);
	}

	MonteCarlo<SquareIsingModel> mc(models);

	auto df = mc.generate_samples(sampler, 10000, 1000, 400, 1, true);
	df.save("data.json");


}
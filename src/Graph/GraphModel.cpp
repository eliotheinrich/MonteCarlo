#include "GraphModel.h"

void GraphModel::init_params(ull N) {
	this->N = N;
}

void GraphModel::init() {
	this->edges = std::vector<std::vector<int>>(N, std::vector<int>(0));
	this->vals = std::vector<int>(N);

	for (int i = 0; i < N; i++) {
		vals[i] = 0;
		for (int j = 0; j < N; j++) {
			edges[i][j] = 0;
		}
	}
}

void GraphModel::generate_mutation() {
	mut.i = rand() % N;
	mut.j = rand() % N;
	while (mut.j == mut.i) {
		mut.j = rand() % N;
	}
}

void GraphModel::accept_mutation() {
	return;
}

void GraphModel::reject_mutation() {
	toggle_edge(mut.i, mut.j);
}

double GraphModel::energy() const {
	double E = 0;

	for (int i = 0; i < N; i++) {
		E += onsite_energy(i);
		E += bond_energy(i);
	}
	return E;
}

double GraphModel::energy_change() {
	double E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
	toggle_edge(mut.i, mut.j);
	double E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

	return E2 - E1;
}

void GraphModel::toggle_edge(int i, int j) {
	edges[i][j] = 1 - edges[i][j];
	edges[j][i] = 1 - edges[j][i];
}

int GraphModel::deg(int i) const {
	int d = 0;
	for (int j = 0; j < N; j++) {
		d += edges[i][j];
	}
	return d;
}

double GraphModel::get_connectivity() const {
	double c = 0.;
	for (int i = 0; i < N; i++) {
		c += deg(i);
	}
	return c/N;
}

std::map<std::string, Sample> GraphModel::take_samples() const {
	std::map<std::string, Sample> samples;
	samples.emplace("connectivity", get_connectivity());
	return samples;
}
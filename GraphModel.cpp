#ifndef GRAPHMC_
#define GRAPHMC_
#include "GraphModel.h"

GraphModel::GraphModel(int N) {
	this->N = N;

	this->edges = std::vector<std::vector<int>>(N, std::vector<int>(0));
	this->vals = std::vector<int>(N);
}

void GraphModel::init() {
	for (int i = 0; i < N; i++) {
		vals[i] = 0;
		for (int j = 0; j < N; j++) {
			edges[i][j] = 0;
		}
	}
}

void GraphModel::generate_mutation() {
	mut.i = std::rand() % N;
	mut.j = std::rand() % N;
	while (mut.j == mut.i) {
		mut.j = std::rand() % N;
	}
}

void GraphModel::accept_mutation() {
	return;
}

void GraphModel::reject_mutation() {
	toggle_edge(mut.i, mut.j);
}

const float GraphModel::energy() {
	float E = 0;

	for (int i = 0; i < N; i++) {
		E += onsite_energy(i);
		E += bond_energy(i);
	}
	return E;
}

const float GraphModel::energy_change() {
	float E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
	toggle_edge(mut.i, mut.j);
	float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

	return E2 - E1;
}

void GraphModel::toggle_edge(int i, int j) {
	edges[i][j] = 1 - edges[i][j];
	edges[j][i] = 1 - edges[j][i];
}

int GraphModel::deg(int i) {
	int d = 0;
	for (int j = 0; j < N; j++) {
		d += edges[i][j];
	}
	return d;
}

float GraphModel::get_connectivity() {
	float c = 0.;
	for (int i = 0; i < N; i++) {
		c += deg(i);
	}
	return c/N;
}

#endif

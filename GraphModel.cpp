#ifndef SPINMC_
#define SPINMC_

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <random>
#include "MonteCarlo.cpp"
#include "Utility.cpp"

class GraphModel : virtual public MCModel {
	// Generic graph model
    private:
        struct GraphMutation {
            int i;
			int j;
        };


    public:
        int N;

        float acceptance;
		std::vector<std::vector<int>> edges;
		std::vector<int> vals;

        // Mutation being considered is stored as an attribute of the model
        GraphMutation mut;

        GraphModel() {}

        GraphModel(int N) {
            this->N = N;

            this->edges = std::vector<int>(std::vector<int>(0, N), N);
			this->vals = std::vector<int>(N);
        }

		void init() {
			for (int i = 0; i < N; i++) {
				vals[i] = 0;
				for (int j = 0; j < N; j++) {
					edges[i][j] = 0;
				}
			}
		}

        void generate_mutation() {
			mut.i = std::rand() % N;
			mut.j = std::rand() % N;
			while (mut.j == mut.i) {
				mut.j = std::rand() % N;
			}
        }

        void accept_mutation() {
			return;
        }

        void reject_mutation() {
			toggle_edge(i, j);
        }

        virtual const float onsite_energy(int i)=0;

        virtual const float bond_energy(int i)=0;

        const float energy() {
            float E = 0;

            for (int i = 0; i < V; i++) {
                E += onsite_energy(i);
                E += bond_energy(i);
            }
            return E;
        }

        const float energy_change() {
            float E1 = onsite_energy(mut.i) + 2*bond_energy(mut.i);
			toggle_edge(mut.i, mut.j);
            float E2 = onsite_energy(mut.i) + 2*bond_energy(mut.i);

            return E2 - E1;
        }

		void toggle_edge(int i, int j) {
			edges[i][j] = 1 - edges[i][j];
			edges[j][i] = 1 - edges[j][i];
		}

		int deg(int i) {
			int d = 0;
			for (int j = 0; j < N; j++) {
				d += edges[i][j];
			}
			return d;
		}

		float get_connectivity() {
			float c = 0.;
			for (int i = 0; i < N; i++) {
				c += deg(i);
			}
			return c/N;
		}
};

#endif

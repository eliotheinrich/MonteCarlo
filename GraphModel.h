#ifndef GRAPHMC_
#define GRAPHMC_

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
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

		GraphModel();
        GraphModel(int N);
		void init();

        void generate_mutation();
        void accept_mutation();
        void reject_mutation();

        virtual const float onsite_energy(int i)=0;
        virtual const float bond_energy(int i)=0;
        const float energy();

        const float energy_change();

		void toggle_edge(int i, int j);
		int deg(int i);

		float get_connectivity();
};

#endif
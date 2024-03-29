#ifndef GRAPH_MC_H
#define GRAPH_MC_H

#include <vector>
#include "MonteCarlo.h"

class GraphModel : virtual public MCModel {
	// Generic graph model
    private:
        struct GraphMutation {
            int i;
			int j;
        };


    public:
        ull N;

        double acceptance;
		std::vector<std::vector<int>> edges;
		std::vector<int> vals;

        // Mutation being considered is stored as an attribute of the model
        GraphMutation mut;

		GraphModel() {}
        virtual ~GraphModel() {}

        void init_params(ull N);
		virtual void init();
        virtual ull system_size() const {
            return N;
        }

        virtual void generate_mutation();
        virtual void accept_mutation();
        virtual void reject_mutation();

        virtual double onsite_energy(int i) const=0;
        virtual double bond_energy(int i) const=0;

        virtual double energy() const;
        virtual double energy_change();

		void toggle_edge(int i, int j);
		int deg(int i) const;

		double get_connectivity() const;
        virtual data_t take_samples() const;
};

#endif
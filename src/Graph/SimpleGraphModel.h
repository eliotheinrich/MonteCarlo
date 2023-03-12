#ifndef SIMPLEGRAPH_H
#define SIMPLEGRAPH_H

#include "GraphModel.h"

class SimpleGraphModel : public GraphModel {
    private:
        int N;
        float J;
        
    public:
        SimpleGraphModel(Params &params);

        virtual MCModel* clone(Params &params) { return new SimpleGraphModel(params); }

        virtual double onsite_energy(int i) const;
        virtual double bond_energy(int i) const;

};

#endif
#ifndef XXZ_H
#define XXZ_H

#include <vector>
#include <Eigen/Dense>
#include "Spin3DModel.h"

#define DEFAULT_LAYERS 1

class XXZHeis : public Spin3DModel {
    private: 
        int N;
        int L;
        float J;
        float K;
        float A;

        bool sample_helicity;

    public:
        XXZHeis(Params &params);

        std::vector<double> vorticity() const;

        virtual double onsite_func(const Eigen::Vector3d &S) const override;

        CLONE(MCModel, XXZHeis)
};

#endif
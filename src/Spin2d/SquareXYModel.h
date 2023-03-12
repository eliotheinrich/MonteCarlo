#ifndef SQUAREXY_H
#define SQUAREXY_H

#include <vector>
#include <Eigen/Dense>
#include "Spin2DModel.h"

#define DEFAULT_LAYERS 1

class SquareXYModel : public Spin2DModel {
    private:
        int N;
        int L;
        int mut_mode;
        float J;
        float B;
        float Bp;
        float Bx;
        float By;

    public:
        SquareXYModel(Params &params);

        virtual MCModel* clone(Params &params) { return new SquareXYModel(params); };

        inline std::vector<double> vorticity() const;

        float p(int i) const;
        float e1() const;
        float e2() const;
        float U2() const;
        std::vector<double> twist_stiffness() const;

        virtual double onsite_func(const Eigen::Vector2d &S) const;

        void over_relaxation_mutation();

        void generate_mutation();
};

#endif
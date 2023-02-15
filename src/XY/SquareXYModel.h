#ifndef SQUAREXY_H
#define SQUAREXY_H

#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include <complex>
#include "../XYModel.h"
#include "../Utility.cpp"


class SquareXYModel : public XYModel {
    public:
        int N;
        int L;
        int mut_mode;
        float J;
        float B;
        float Bp;
        float Bx;
        float By;

        SquareXYModel(int N, int L, float J, float B, float Bp);

        SquareXYModel* clone();

        inline std::vector<double> vorticity();

        float p(int i);
        float e1();
        float e2();
        float U2();
        std::vector<double> twist_stiffness();

        virtual double onsite_func(const Eigen::Vector2d &S) const;

        void over_relaxation_mutation();

        void generate_mutation();
};

#endif
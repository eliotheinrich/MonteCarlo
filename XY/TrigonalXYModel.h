#ifndef TRIGONALXY_H
#define TRIGONALXY_H

#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../XYModel.h"
#include "../Utility.cpp"

class TrigonalXYModel : public XYModel {
    public:
        int N;
        int L;
        float J;
        float A;

        int mut_mode;

        TrigonalXYModel(int N, int L, float J, float A);

        TrigonalXYModel* clone();

        inline std::vector<double> vorticity();

        virtual double onsite_func(const Eigen::Vector2d &S) const;

        void over_relaxation_mutation();
        void generate_mutation();
};


#endif
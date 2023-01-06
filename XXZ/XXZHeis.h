#ifndef XXZ_H
#define XXZ_H

#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include "../SpinModel.h"
#include "../Utility.cpp"

class XXZHeis : public SpinModel {
    public:
        int N;
        int L;
        float J;
        float K;

        XXZHeis(int N, int L, float J, float K);

        XXZHeis* clone();

        inline std::vector<double> vorticity();

        virtual double onsite_func(const Eigen::Vector3d &S) const;
};

#endif
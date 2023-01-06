#ifndef TRIGONALMODEL_H
#define TRIGONALMODEL_H

#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <Eigen/Dense>
#include <set>
#include "../SpinModel.h"
#include "../MonteCarlo.h"
#include "../Utility.cpp"

class TrigonalModel : public SpinModel {
    public:
        int N;
        int L;
        float J1;
        float J2;
        float K1;
        float K2;
        float K3;
        Eigen::Vector3d B; 

        Eigen::Matrix3d R;
        int mut_counter;
        int mut_mode;

        GaussianDist *dist_r;


    TrigonalModel(int N, int L, float J1, float J2, float K1, float K2, float K3,
                                Eigen::Vector3d B);

    TrigonalModel* clone();

    void over_relaxation_mutation();

    void generate_mutation();

    double onsite_func(const Eigen::Vector3d &S) const;

    // --- FOR DYNAMIC UPDATES --- //
    void rotate_spin(int i, Eigen::Vector3d v, float p);

    const Eigen::Vector3d molecular_field(int i);

    void dynamic_step(float dt);

    // ------ //

    Eigen::Vector3d order_param();

    std::vector<double> twist_derivatives();

    inline double winding(double p1, double p2, double p3);
    std::vector<std::vector<double>> find_vortices();
    std::vector<double> vorticity_correlation(int i);
    std::vector<double> vortex_data();
    double rel_pos(int i, int j);
    std::vector<double> vortex_pairings();
};

#endif
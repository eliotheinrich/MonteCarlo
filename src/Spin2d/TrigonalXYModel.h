#ifndef TRIGONALXY_H
#define TRIGONALXY_H

#include <vector>
#include <Eigen/Dense>
#include "Spin2DModel.h"
    
#define DEFAULT_LAYERS 1

class TrigonalXYModel : public Spin2DModel {
    private:
        int N;
        int L;
        float J;
        float A;

        int mut_mode;

    public:

        TrigonalXYModel(Params &params);

        inline std::vector<double> vorticity() const;

        virtual double onsite_func(const Eigen::Vector2d &S) const override;

        void over_relaxation_mutation();
        void generate_mutation();

        CLONE(MCModel, TrigonalXYModel)
};


#endif
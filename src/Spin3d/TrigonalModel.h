#ifndef TRIGONALMODEL_H
#define TRIGONALMODEL_H

#include <vector>
#include "Spin3DModel.h"
#include "MonteCarlo.h"

#define DEFAULT_LAYERS 1

#define DEFAULT_SAMPLE_INTENSITY true

#define DEFAULT_SAMPLE_LAYER_MAGNETIZATION false

class TrigonalModel : public Spin3DModel {
    private:
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

        bool sample_layer_magnetization;
        bool sample_intensity;

    public:
        TrigonalModel(Params &params);

        void over_relaxation_mutation();

        void generate_mutation();

        double onsite_func(const Eigen::Vector3d &S) const;

        // --- FOR DYNAMIC UPDATES --- //
        void rotate_spin(int i, Eigen::Vector3d v, float p);

        Eigen::Vector3d molecular_field(int i) const;

        void dynamic_step(float dt);
        // ------ //

        // For computing structure factor
        Eigen::Vector3d rel_pos(uint i) const;
        double intensity(Eigen::Vector3d Q) const;

        void add_intensity_samples(std::map<std::string, Sample> &samples) const;
        void add_layer_magnetization_samples(std::map<std::string, Sample> &samples) const;

        virtual std::map<std::string, Sample> take_samples();

        CLONE(MCModel, TrigonalModel)
};

#endif
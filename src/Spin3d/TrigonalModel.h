#ifndef TRIGONALMODEL_H
#define TRIGONALMODEL_H

#include <vector>
#include "Spin3DModel.h"
#include "MonteCarlo.h"

#define DEFAULT_LAYERS 1

#define DEFAULT_SAMPLE_LAYER_MAGNETIZATION false

#define DEFAULT_SAMPLE_INTENSITY false
#define DEFAULT_MAX_L 1.
#define DEFAULT_MIN_L 0.
#define DEFAULT_INTENSITY_RESOLUTION 30


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

        bool sample_magnetization;

        bool sample_helicity;

        bool sample_layer_magnetization;
        bool sample_intensity;
        float max_L;
        float min_L;
        uint intensity_resolution;

    public:
        TrigonalModel(Params &params);

        void over_relaxation_mutation();
        virtual void generate_mutation() override;

        virtual double onsite_func(const Eigen::Vector3d &S) const override;

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

        virtual std::map<std::string, Sample> take_samples() override;

        CLONE(MCModel, TrigonalModel)
};

#endif
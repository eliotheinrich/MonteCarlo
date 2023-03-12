#ifndef SPIN2D_MC_H
#define SPIN2D_MC_H

#include <functional>
#include <vector>
#include <unordered_set>
#include <stack>
#include <Eigen/Dense>
#include "MonteCarlo.h"

class Spin2DModel : virtual public MCModel {
    // Generic 2d spin model in up to 3d lattice
    private:
        // A mutation consists of a change in spin dS on site (n1,n2,n3,s)
        // dS must conserve the norm of S[n1,n2,n3,s]
        struct Spin2DMutation {
            int i;
            Eigen::Vector2d dS;
        };

        struct Spin2DBond {
            int d1;
            int d2;
            int d3;
            int ds;
            Eigen::Vector3d v;
            std::function<float(Eigen::Vector2d, Eigen::Vector2d)> bondfunc;
        };



    public:
        int sl;
        int N1;
        int N2;
        int N3;
        int V;

        float acceptance;
        float sigma;
        std::vector<Eigen::Vector2d> spins;
        std::vector<std::vector<int>> neighbors;
        std::vector<Spin2DBond> bonds;

        static constexpr double alpha = 0.01;
        std::vector<Eigen::Matrix2d> R1s;
        std::vector<Eigen::Matrix2d> R2s;
        std::vector<Eigen::Matrix2d> R3s;


        std::unordered_set<int> s;
        Eigen::Matrix2d s0;

        // Mutation being considered is stored as an attribute of the model
        Spin2DMutation mut;

        Spin2DModel() {}

        void init_params(int sl, int N1, int N2, int N3);

        virtual void init();

        inline int flat_idx(int n1, int n2, int n3, int s) const {
            return n1 + N1*(n2 + N2*(n3 + N3*s));
        }

        inline Eigen::Vector4i tensor_idx(int i) const {
            int n1 = i % N1;
            i = i / N1;
            int n2 = i % N2;
            i = i / N2;
            int n3 = i % N3;
            i = i / N3;
            int s = i % sl;
            Eigen::Vector4i v; v << n1, n2, n3, s;
            return v;
        }

        void randomize_spins();

        void add_bond(int d1, int d2, int d3, int ds, Eigen::Vector3d v, std::function<float(Eigen::Vector2d, Eigen::Vector2d)> bondfunc);

        std::vector<double> twist_stiffness() const;

        Eigen::Vector2d get_magnetization() const;
        
        virtual double onsite_energy(int i) const;
        virtual double bond_energy(int i) const;

        void metropolis_mutation();
        void cluster_update();

        virtual void generate_mutation();
        virtual void accept_mutation();
        virtual void reject_mutation();

        virtual double energy() const;
        virtual double energy_change();

        virtual double onsite_func(const Eigen::Vector2d& S) const = 0;
        
        virtual std::map<std::string, Sample> take_samples() const;

        void save_spins(std::string filename);
};

#endif

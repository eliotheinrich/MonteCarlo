#include "GraphModel.h"

class SimpleGraphModel : public GraphModel {
  public:
    SimpleGraphModel(dataframe::ExperimentParams &params, uint32_t num_threads);

    virtual double onsite_energy(uint32_t i) const override;
    virtual double bond_energy(uint32_t i) const override;

  private:
    uint32_t N;
    double J;
};

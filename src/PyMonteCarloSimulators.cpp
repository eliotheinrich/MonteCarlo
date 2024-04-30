#include "SimpleGraphModel.h"

#include "SquareIsingModel.h"

#include "LDPCIsingModel.h"


#include "SquareXYModel.h"
#include "TrigonalXYModel.h"

#include "XXZHeis.h"
#include "TrigonalModel.h"

#include <PyDataFrame.hpp>

NB_MODULE(montecarlo_bindings, m) {
  // Graph
  EXPORT_SIMULATOR_DRIVER(SimpleGraphModel);

  // Ising
  EXPORT_SIMULATOR_DRIVER(SquareIsingModel);
  EXPORT_SIMULATOR_DRIVER(LDPCIsingModel);

  // Clock
  //EXPORT_SIMULATOR_DRIVER(SquareClockModel<4>);

  // Spin2D
  EXPORT_SIMULATOR_DRIVER(SquareXYModel);
  EXPORT_SIMULATOR_DRIVER(TrigonalXYModel);

  // Spin3D
  EXPORT_SIMULATOR_DRIVER(XXZHeis);
  EXPORT_SIMULATOR_DRIVER(TrigonalModel);
}



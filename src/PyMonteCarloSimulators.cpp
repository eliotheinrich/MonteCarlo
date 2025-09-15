#include "Graph/SimpleGraphModel.h"

#include "Ising/SquareIsingModel.h"

#include "MultibodyIsing/LDPCIsingModel.h"

#include "Spin2d/SquareXYModel.h"
#include "Spin2d/TrigonalXYModel.h"

#include "Spin3d/XXZHeis.h"
#include "Spin3d/TrigonalModel.h"
#include "Spin3d/AltermagnetModel.h"
#include "Spin3d/HelixModel.h"

#include <PyDataFrame.hpp>
#include <PyUtils.hpp>

NB_MODULE(montecarlo_bindings, m) {
  // Graph
  EXPORT_SIMULATOR(SimpleGraphModel);

  // Ising
  EXPORT_SIMULATOR(SquareIsingModel);
  EXPORT_SIMULATOR(LDPCIsingModel);

  // Clock
  //EXPORT_SIMULATOR(SquareClockModel<4>);

  // Spin2D
  EXPORT_SIMULATOR(SquareXYModel);
  EXPORT_SIMULATOR(TrigonalXYModel);

  // Spin3D
  EXPORT_SIMULATOR(XXZHeis);
  EXPORT_SIMULATOR(TrigonalModel);
  EXPORT_SIMULATOR(AltermagnetModel)
    .def("to_graph", [](AltermagnetModel& self) { 
      auto graph = self.to_graph();
      return std::make_tuple(graph.edges, graph.vals);
    });

  EXPORT_SIMULATOR(HelixModel);
}



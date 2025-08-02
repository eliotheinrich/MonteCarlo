#include "SimpleGraphModel.h"

#include "SquareIsingModel.h"

#include "LDPCIsingModel.h"

#include "SquareXYModel.h"
#include "TrigonalXYModel.h"

#include "XXZHeis.h"
#include "TrigonalModel.h"
#include "Selenium.h"

#include <PyDataFrame.hpp>
#include <PyQutils.hpp>

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
  EXPORT_SIMULATOR(SeleniumModel)
    .def("to_graph", [](SeleniumModel& self) { 
      auto graph = self.to_graph();
      return std::make_tuple(graph.edges, graph.vals);
    });
}



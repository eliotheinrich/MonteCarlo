add_library(mc_graph
	GraphModel.cpp
	SimpleGraphModel.cpp
)
target_include_directories(mc_graph PUBLIC ${CMAKE_SOURCE_DIR}/src)
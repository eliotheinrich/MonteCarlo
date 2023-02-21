add_library(graph
	GraphModel.cpp GraphModel.h
	SimpleGraphModel.cpp SimpleGraphModel.h
)
target_link_libraries(graph montecarlo)
target_include_directories(graph PRIVATE ${CMAKE_SOURCE_DIR}/src $ENV{EIGEN_DIR})
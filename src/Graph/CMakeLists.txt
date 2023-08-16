add_library(graph
	GraphModel.cpp
	SimpleGraphModel.cpp
)
target_link_libraries(graph montecarlo)
target_include_directories(graph PUBLIC ${CMAKE_SOURCE_DIR}/src)
add_library(ising
	IsingModel.cpp
	SquareIsingModel.cpp
)
target_link_libraries(ising montecarlo)
target_include_directories(ising PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_compile_definitions(ising PUBLIC CLUSTER_UPDATE)

add_library(ising_metropolis
	IsingModel.cpp
	SquareIsingModel.cpp
)
target_link_libraries(ising_metropolis montecarlo)
target_include_directories(ising_metropolis PUBLIC ${CMAKE_SOURCE_DIR}/src)

include_directories($ENV{EIGEN_DIR} $ENV{CTPL_DIR} ${Boost_INCLUDE_DIR})



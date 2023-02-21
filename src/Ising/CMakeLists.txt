add_library(ising
	IsingModel.cpp IsingModel.h
	SquareIsingModel.cpp SquareIsingModel.h
)
target_link_libraries(ising montecarlo)
target_include_directories(ising PRIVATE ${CMAKE_SOURCE_DIR}/src $ENV{EIGEN_DIR})

add_executable(run_ising RunIsing.cpp)
target_link_libraries(run_ising LINK_PUBLIC ising montecarlo)
target_include_directories(run_ising PRIVATE ${CMAKE_SOURCE_DIR}/src $ENV{EIGEN_DIR})


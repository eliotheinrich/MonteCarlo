add_library(ising
	IsingModel.cpp
	SquareIsingModel.cpp
)
target_link_libraries(ising montecarlo)
target_include_directories(ising PUBLIC ${CMAKE_SOURCE_DIR}/src ${EIGEN_DIR} ${DF_DEPENDENCIES})

include_directories($ENV{EIGEN_DIR} $ENV{CTPL_DIR} ${Boost_INCLUDE_DIR})
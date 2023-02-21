add_library(heisenberg
	SpinModel.cpp SpinModel.h
	Trigonal/TrigonalModel.cpp Trigonal/TrigonalModel.h
	XXZ/XXZHeis.cpp XXZ/XXZHeis.h
)
target_link_libraries(heisenberg montecarlo)
target_include_directories(heisenberg PRIVATE ${CMAKE_SOURCE_DIR}/src $ENV{EIGEN_DIR})
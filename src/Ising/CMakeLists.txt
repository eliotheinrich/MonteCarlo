add_library(ising
	IsingModel.cpp
	SquareIsingModel.cpp
)
target_include_directories(ising PUBLIC ${CMAKE_SOURCE_DIR}/src)

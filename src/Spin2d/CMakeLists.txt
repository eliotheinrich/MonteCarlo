add_library(spin2d
	Spin2DModel.cpp
	SquareXYModel.cpp
	TrigonalXYModel.cpp
)
target_link_libraries(spin2d montecarlo)
target_include_directories(spin2d PUBLIC ${CMAKE_SOURCE_DIR}/src)
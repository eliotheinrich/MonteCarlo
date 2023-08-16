add_library(spin3d
	Spin3DModel.cpp
	TrigonalModel.cpp
	XXZHeis.cpp
)
target_link_libraries(spin3d montecarlo)
target_include_directories(spin3d PUBLIC ${CMAKE_SOURCE_DIR}/src)
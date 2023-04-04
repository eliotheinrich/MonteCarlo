add_library(spin2d
	Spin2DModel.cpp
	SquareXYModel.cpp
	TrigonalXYModel.cpp
)
target_link_libraries(spin2d montecarlo)
target_include_directories(spin2d PUBLIC ${CMAKE_SOURCE_DIR}/src ${EIGEN_DIR} ${DF_DEPENDENCIES})
target_compile_definitions(spin2d PUBLIC CLUSTER_UPDATE)

add_library(spin2d_metropolis
	Spin2DModel.cpp
	SquareXYModel.cpp
	TrigonalXYModel.cpp
)
target_link_libraries(spin2d_metropolis montecarlo)
target_include_directories(spin2d_metropolis PUBLIC ${CMAKE_SOURCE_DIR}/src ${EIGEN_DIR} ${DF_DEPENDENCIES})
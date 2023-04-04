add_library(spin3d
	Spin3DModel.cpp
	TrigonalModel.cpp
	XXZHeis.cpp
)
target_link_libraries(spin3d montecarlo)
target_include_directories(spin3d PUBLIC ${CMAKE_SOURCE_DIR}/src ${EIGEN_DIR} ${DF_DEPENDENCIES})
target_compile_definitions(spin3d PUBLIC CLUSTER_UPDATE)

add_library(spin3d_metropolis
	Spin3DModel.cpp
	TrigonalModel.cpp
	XXZHeis.cpp
)
target_link_libraries(spin3d_metropolis montecarlo)
target_include_directories(spin3d_metropolis PUBLIC ${CMAKE_SOURCE_DIR}/src ${EIGEN_DIR} ${DF_DEPENDENCIES})
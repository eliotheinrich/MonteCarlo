find_package(glfw3 3.4)

set(GLAD_PATH $ENV{GLAD_PATH})
include_directories("${GLAD_PATH}/include")

set(SOURCE_FILES
  Animator.cpp
  Shader.cpp
)

# Add either glad.c or glad.cpp to the source files list
if(EXISTS "${GLAD_PATH}/src/glad.c")
  list(APPEND SOURCE_FILES "${GLAD_PATH}/src/glad.c")
elseif(EXISTS "${GLAD_PATH}/src/glad.cpp")
  list(APPEND SOURCE_FILES "${GLAD_PATH}/src/glad.cpp")
else()
  message(FATAL_ERROR "Neither glad.c nor glad.cpp was found in ${GLAD_PATH}/src")
endif()

message("Inside GLFW src. SOURCE_FILES = ${SOURCE_FILES}")

add_library(simulator_display ${SOURCE_FILES})
#target_compile_definitions(simulator_display  PRIVATE SIMULATOR_DISPLAY_PATH="${CMAKE_SOURCE_DIR}/src/SimulatorDisplay")
target_link_libraries(simulator_display PRIVATE glfw)

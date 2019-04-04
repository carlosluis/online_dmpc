cmake_minimum_required(VERSION 3.0)

project(online_planning)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
find_package(Eigen3 REQUIRED)
include_directories(
	${PROJECT_SOURCE_DIR}/include
	${EIGEN3_INCLUDE_DIR}
)

message("-- Eigen library directory: ${EIGEN3_INCLUDE_DIR}")
message("-- Source path: ${PROJECT_SOURCE_DIR}")

add_executable(run
	src/main.cpp
)

add_library(bezier
	src/bezier.cpp)

add_library(model
	src/model.cpp)

add_library(generator
	src/generator.cpp)

target_link_libraries(generator
	bezier
	pthread
	eigen-quadprog
	model)

target_link_libraries(run
	generator
)


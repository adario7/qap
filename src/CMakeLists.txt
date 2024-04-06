cmake_minimum_required(VERSION 3.20)

project(qap)

# genereate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# find sources
file(GLOB_RECURSE LIB CONFIGURE_DEPENDS "lib/*.cc" "include/*.hh")
add_library(lib ${LIB})
include_directories(include)
set(CMAKE_CXX_STANDARD 17)

# enable sanitizers in debug mode
option(SANITIZE "Enable sanitizaers in debug builds" ON)
if(SANITIZE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
	add_compile_options(-fomit-frame-pointer -fsanitize=address,undefined)
	add_link_options(-fomit-frame-pointer -fsanitize=address,undefined)
endif()

# dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
find_package(CPLEX REQUIRED)
include_directories(${CPLEX_INCLUDE_DIRS})

# binary targets
add_executable(toff toff.cc)
target_link_libraries(toff PRIVATE lib ${CPLEX_LIBRARIES})

add_executable(cb cb.cc)
target_link_libraries(cb PRIVATE lib ${CPLEX_LIBRARIES})

add_executable(nocb experiments/nocb.cc)
target_link_libraries(nocb PRIVATE lib ${CPLEX_LIBRARIES})

add_executable(custom experiments/custom.cc)
target_link_libraries(custom PRIVATE lib ${CPLEX_LIBRARIES})

add_executable(mixed experiments/mixed.cc)
target_link_libraries(mixed PRIVATE lib ${CPLEX_LIBRARIES})

add_executable(a_star experiments/a_star.cc)
target_link_libraries(a_star PRIVATE lib ${CPLEX_LIBRARIES})

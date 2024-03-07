cmake_minimum_required(VERSION 3.20)

project(qap)

# genereate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# find sources
file(GLOB_RECURSE LIB CONFIGURE_DEPENDS "lib/*.cc" "include/*.hh")
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
add_executable(toff ${LIB} toff.cc)
target_link_libraries(toff PRIVATE ${CPLEX_LIBRARIES})

add_executable(nocb ${LIB} nocb.cc)
target_link_libraries(nocb PRIVATE ${CPLEX_LIBRARIES})

add_executable(cb ${LIB} cb.cc)
target_link_libraries(cb PRIVATE ${CPLEX_LIBRARIES})

add_executable(cb2 ${LIB} experiments/cb2.cc)
target_link_libraries(cb2 PRIVATE ${CPLEX_LIBRARIES})

add_executable(custom ${LIB} experiments/custom.cc)
target_link_libraries(custom PRIVATE ${CPLEX_LIBRARIES})

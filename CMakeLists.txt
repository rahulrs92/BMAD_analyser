cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
#set(CMAKE_C_COMPILER "gcc-14")
#set(CMAKE_CXX_COMPILER "g++-14")

project(bmad_analyser)

find_package(ROOT REQUIRED COMPONENTS Graf Gpad Hist Tree TreePlayer)
include(${ROOT_USE_FILE})

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB lv1 ${PROJECT_SOURCE_DIR}/lv1/*.cc)
file(GLOB lv2 ${PROJECT_SOURCE_DIR}/lv2/*.cc)
file(GLOB analyser ${PROJECT_SOURCE_DIR}/analyser/*.cc)
file(GLOB acquisition ${PROJECT_SOURCE_DIR}/acquisition/*.cc)

add_executable(bmad_analyser ${sources} ${headers} ${lv1} ${lv2} ${analyser} ${acquisition} main.cpp)
install(TARGETS bmad_analyser RUNTIME DESTINATION bin)
target_link_libraries(bmad_analyser ROOT::Core ROOT::Graf ROOT::Gpad ROOT::Hist ROOT::Tree ROOT::TreePlayer)
target_compile_features(bmad_analyser PRIVATE cxx_std_17)
target_include_directories(bmad_analyser PRIVATE external/Headers)

include_directories(${PROJECT_SOURCE_DIR}/include)



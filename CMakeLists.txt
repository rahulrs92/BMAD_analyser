cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(bmad_analyser)

find_package(ROOT REQUIRED COMPONENTS Graf Gpad Hist)
include(${ROOT_USE_FILE})

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB lv1 ${PROJECT_SOURCE_DIR}/lv1/*.cc)
file(GLOB lv2 ${PROJECT_SOURCE_DIR}/lv2/*.cc)
file(GLOB analyser ${PROJECT_SOURCE_DIR}/analyser/*.cc)

add_executable(bmad_analyser ${sources} ${headers} ${lv1} ${lv2} ${analyser} main.cpp)
install(TARGETS bmad_analyser RUNTIME DESTINATION bin)
target_link_libraries(bmad_analyser ROOT::Core ROOT::Graf ROOT::Gpad ROOT::Hist)
target_compile_features(bmad_analyser PRIVATE cxx_std_17)

include_directories(${PROJECT_SOURCE_DIR}/include)



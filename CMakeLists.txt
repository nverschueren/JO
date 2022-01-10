project(jos)
set(project_SOURCES jos.cpp)
set(CMAKE_CXX_FLAGS "-Ofast -Wall -Wextra")
add_executable(jos ${project_SOURCES})
target_link_libraries(jos armadillo glut GL GLU)


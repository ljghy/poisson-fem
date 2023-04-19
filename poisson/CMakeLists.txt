add_executable(poisson)

target_sources(poisson
PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Cell.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver.cpp

)

target_include_directories(poisson
PRIVATE 
)

find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(poisson
PRIVATE 
    Eigen3::Eigen
)

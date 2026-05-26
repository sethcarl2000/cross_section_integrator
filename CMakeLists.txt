#---------------------------------------------------------
#   I copied most of this CMakeLists.txt from: https://cliutils.gitlab.io/modern-cmake/examples/root-dict/README.html 
cmake_minimum_required(VERSION 3.11)

#---------------------------------------------------------
#   Define the project 
project(cs_integrator VERSION 0.1
        LANGUAGES CXX)  

#---------------------------------------------------------
#   Define this preprocessor macro to enable debug statements 
option(DEBUG "Enable debug printing" OFF)
if (DEBUG)
    add_compile_definitions(DEBUG)
endif()

#---------------------------------------------------------
#   Define this preprocessor macro to enable (more expensive) range-checks 


#---------------------------------------------------------
#   Which include directories to use. 
set(CMAKE_CXX_STANDARD
    20
    CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

#---------------------------------------------------------
#   Find ROOT (necessary)
find_package(ROOT 6.20 REQUIRED COMPONENTS Core RIO Physics)
find_package(GSL REQUIRED)

#---------------------------------------------------------
#   Define sources and headers
set(sources 
    processes/trident.cxx
    processes/bh_photoproduction.cxx
    processes/elastic.cxx
    PolarFourVec.cxx
    MetropolisIntegrator.cxx
    M2_expr.cxx
    DifferentialCSIntegrator.cxx

)

set(headers
    IntegSettings.hxx
    Bound.hxx
    PolarFourVec.hxx
    FourVec.hxx
    processes.hxx 
    DifferentialCSIntegrator.hxx
    IntegrationResult.hxx
    some_numbers.hxx
    MetropolisIntegrator.hxx
    argparse.hpp
)

#---------------------------------------------------------
#   Add application executable
#add_executable(main main.cxx ${sources} ${headers})
#target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
#target_link_libraries(main PRIVATE ${ROOT_LIBRARIES})

#add_executable(avg_dist avg_dist.cxx ${sources} ${headers})
#target_include_directories(avg_dist PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
#target_link_libraries(avg_dist PRIVATE ${ROOT_LIBRARIES})

#---------------------------------------------------------
#   Add Libraries
add_executable(cs_integrator cs_integrator.cxx ${sources} ${headers})
target_include_directories(cs_integrator PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(cs_integrator PRIVATE -g)
target_link_libraries(cs_integrator PUBLIC ROOT::Core Threads::Threads)
target_compile_options(cs_integrator PRIVATE -pthread)

# add threads & other args 
find_package(Threads REQUIRED)
target_link_libraries(cs_integrator PUBLIC Threads::Threads)

target_include_directories(cs_integrator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})
target_link_libraries(cs_integrator PRIVATE ${ROOT_LIBRARIES} ${GSL_LIBRARIES})


#---------------------------------------------------------
#   Generate the ROOT dictionary
#root_generate_dictionary(G__cs_integrator ${headers} LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries
#add_library(cs_integrator SHARED ${headers} ${sources} G__cs_integrator.cxx)
#target_include_directories(cs_integrator PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
#target_link_libraries(cs_integrator PUBLIC ROOT::Core)


#target_include_directories(pick_sieveholes PRIVATE
#    "${CMAKE_CURRENT_SOURCE_DIR}"
#    "${ROOT_INCLUDE_DIRS}"
#)

#target_link_libraries(pick_sieveholes PRIVATE
#    ${ROOT_LIBRARIES}
#    cs_integrator
#)

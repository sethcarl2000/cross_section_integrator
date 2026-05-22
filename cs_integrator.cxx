#include "FourVec.hxx"
#include "PolarFourVec.hxx"
#include "processes.hxx"
#include "argparse.hpp"

#include <iostream>

//#define EDCS_DEBUG
#include "EstimateDifferentialCS.hxx"

namespace {
    constexpr int D = 4; 

    constexpr double deg = 3.1415926536 / 180.; 

    constexpr double me = 0.501;
}

int main(int argc, char* argv[])
{
    std::cout << "attempting integration..." << std::endl; 

    const double beam_E = 2205.;

    PolarFourVec P0, P1;
    
    P0.cos_theta = 1.0;
    P0.phi       = 0.; 
    P0.energy    = beam_E; 
    P0.mass2     = me*me; 
    
    P1.cos_theta = 1.0;
    P1.phi       = 0.; 
    P1.energy    = beam_E/2.; 
    P1.mass2     = me*me; 
    
    double amp = EstimateDifferentialCS<4>(
        processes::trident,
        P0, 
        P1, 1, { me, me },
        1e6 
    );

    std::printf("done. final amplitude: %e\n", amp); 
    
    return 0; 
}
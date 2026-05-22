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
    
    P1.cos_theta = 0.95;
    P1.phi       = 0.; 
    P1.energy    = beam_E*0.95; 
    P1.mass2     = me*me; 
    
    int n_samples = 20;

    double amp{0.}, amp2{0.};
    for (int i=0; i<n_samples; i++) {
        
        double amp_i = EstimateDifferentialCS<3>(
            processes::bh_photoproduction,
            P0, 
            P1, 1, { 0. },
            3.e7, 
            Setting::kVerbose, Setting::kMISER
        );

        amp  += amp_i;
        amp2 += amp_i*amp_i; 
    }   
    
    double variance = amp2/((double)n_samples) - std::pow( amp/((double)n_samples), 2 );

    std::printf("done. final amplitude: %.5e  +/-  %.4e\n", amp/((double)n_samples), std::sqrt(variance)); 
    
    return 0; 
}
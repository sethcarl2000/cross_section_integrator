#include <iostream>
#include <stdio.h>
#include "CrossSectionIntegrator.hxx"
#include "FourVec.hxx"
#include <functional> 
#include <array> 
#include "processes.hxx"

namespace {

    constexpr int D = 4; 

    //mass of electron MeV. 
    constexpr double me = 0.501; 

    //norm2 of the spacelike part of a 3-vector
    double space_norm2(const FourVec& v) {
        return v.get(1)*v.get(1) + v.get(2)*v.get(2) + v.get(3)*v.get(3); 
    }

    ///
    /// P[0] - incoming momenta
    /// P[1] - outgoing momenta 
    double Rutherford_Scattering(const std::array<FourVec,2>& P) {

        double Q2 = space_norm2(P[1] - P[0]); 

        double tr = 4.*P[0].norm2() + 8.*P[0].get(0)*P[1].get(0) - 4.*(P[0]*P[1]); 
        return tr / ( Q2*Q2 + 1. );  
    }
}

int test()
{
    std::cout << "attempting integration..." << std::endl; 

    const double beam_E = 2200.;

    FourVec P0{ std::sqrt(beam_E*beam_E + me*me), 0., 0., +beam_E }; 

    FourVec e_rest{ me, 0., 0., 0. }; 

    FourVec P1 = Nudge( e_rest, 0,  20,  1998. );
    FourVec Pp = Nudge( e_rest, 0., 70.7, 70.7 );
    FourVec Pm = Nudge( e_rest, 0., -70.7, 70.7 );

    std::array<FourVecParameters,4> inputs_trident = {
        FourVecParameters{ .is_fixed=true,   .val=P0, .momenta_scan_amplitude=beam_E/2500., .mass=me },
        FourVecParameters{ .is_fixed=true,   .val=P1, .momenta_scan_amplitude=beam_E/1000., .mass=me },
        FourVecParameters{ .is_fixed=false,  .val=Pp, .momenta_scan_amplitude=beam_E/2500., .mass=me },
        FourVecParameters{ .is_fixed=false,  .val=Pm, .momenta_scan_amplitude=beam_E/2500., .mass=me }
        //,
    //    FourVecParameters{ .is_fixed=false, .val=P0, .momenta_scan_amplitude=1., .mass=me }
    };

    CrossSectionIntegrator<D>(processes::trident, inputs_trident, 1e6, 1e4, Setting::kConserveEnergy ); 

    std::cout << "done." << std::endl; 
    return 0; 
}   
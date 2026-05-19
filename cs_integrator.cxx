#include "FourVec.hxx"
#include "CrossSectionIntegrator.hxx"
#include "processes.hxx"
#include "argparse.hpp"

#include <iostream>

namespace {
    constexpr int D = 4; 

    constexpr double deg = 3.1415926536 / 180.; 

    constexpr double me = 0.501;
}

int main(int argc, char* argv[])
{
    std::cout << "attempting integration..." << std::endl; 

    const double beam_E = 2200.;

    FourVec P0{ std::sqrt(beam_E*beam_E + me*me), 0., 0., +beam_E }; 

    FourVec e_rest{ me, 0., 0., 0. }; 

    FourVec P1 = Nudge( e_rest, 0,  2000.*std::sin(5.*deg), 2000.*std::cos(5.*deg)  );
    FourVec Pp = Nudge( e_rest, 0., +100*std::sin(15*deg),    100*std::cos(15*deg)    );
    FourVec Pm = Nudge( e_rest, 0., -100*std::sin(15*deg),    100*std::cos(15*deg)    );

    std::array<FourVecParameters,4> inputs_trident = {
        FourVecParameters{ .is_fixed=true,   .val=P0, .momenta_scan_amplitude=0., .mass=me },
        FourVecParameters{ .is_fixed=true,   .val=P1, .momenta_scan_amplitude=0., .mass=me },
        FourVecParameters{ .is_fixed=false,  .val=Pp, .momenta_scan_amplitude=1e-2, .mass=me },
        FourVecParameters{ .is_fixed=false,  .val=Pm, .momenta_scan_amplitude=1e-2, .mass=me }
        //,
    //    FourVecParameters{ .is_fixed=false, .val=P0, .momenta_scan_amplitude=1., .mass=me }
    };

    long double amp = CrossSectionIntegrator<D>(
        processes::trident, 
        inputs_trident, 
        10e6, 
        1e5, 
        Setting::kConserveEnergy | Setting::kAutoAdjustScan, 
        2.5e6
    ); 

    std::printf("done. final amplitude: %Lf\n", amp); 
    
    return 0; 
}
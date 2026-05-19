#ifndef CrossSectionIntegrator_H
#define CrossSectionIntegrator_H

#include <cmath> 
#include <stdio.h>
#include <functional> 
#include <random> 
#include <array> 
#include <iostream> 
#include <limits> 
#include "FourVec.hxx"

struct FourVecParameters {
    bool is_fixed; //if 'true', then it will not be integrated over
    FourVec val; //starting value 
    double momenta_scan_amplitude; //amplitude of the momentum to scan
    double mass; //mass of the particle (MeV/c^2)
};

namespace Setting {
    enum Bit : int {
        kNone = 0,
        kConserveEnergy = 1 << 0, //enforce external energy conservation for each update
        kConserveMomenta = 1 << 1, //enforce external momentum conservation for each update 
        kAutoAdjustScan = 1 << 2, //automatically adjust the scan rate to keep the acceptance prob around ~50% 
    };
}

//integrates a square matrix element over external momenta
template<int D> double CrossSectionIntegrator(
    std::function<double(const std::array<FourVec,D>&)> fcn, 
    const std::array<FourVecParameters,D> momenta_inputs, 
    long long int n_steps, 
    long long int steps_between_reports,
    Setting::Bit mode = Setting::kConserveEnergy
)
{       
    const double kNaN = std::numeric_limits<double>::quiet_NaN(); 

    /// @return 'true' if arg is nan
    auto is_nan = [](double x) { return x!=x; };

    //initialize random-number generators
    std::mt19937 twister; 
    std::uniform_real_distribution<double> runiform{0., 1.}; 
    auto Rndm = std::bind(runiform, twister); 

    std::normal_distribution<double> rgaus{0., 1.};
    auto Gaus = std::bind(rgaus, twister); 

    //the (relativistic) integral measure for 4-vectors (which are on-shell)
    std::array<std::function<double(FourVec)>,D> measure;

    //the momenta we will vary
    std::array<FourVec,D> momenta; 

    //momenta to integrate (rest are fixed)
    std::vector<int> momenta_to_integrate; 

    for (int i=0; i<D; i++) {

        const auto& inp = momenta_inputs[i];

        momenta[i] = inp.val;
        if (!inp.is_fixed) { momenta_to_integrate.push_back(i); }

        double m2 = inp.mass*inp.mass; 
        measure[i] = [m2](const FourVec& p){ return 1./std::sqrt( m2 + p.norm2() ); }; 
    }
    std::cout << " total momenta to integrate over: " << momenta_to_integrate.size() << std::endl; 

    //try to perform monte-carlo integration of variables    
    long double amp_sum = 0.; 

    //get the starting value of the amplitude 
    double amp = fcn(momenta);

    if (is_nan(amp) || amp < 0.) {
        std::cerr << "<"<<__func__<<">: Fatal error: ~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "   amplitude is negative / NaN. amp=" << amp << "\n"
        "   external momenta: \n"; 
        for (int i=0; i<D; i++) {
            std::cerr << "["<<i<<"]: " << momenta[i] << "\n"; 
        }
        std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        return -kNaN; 
    }

    //if all momenta are fixed, then just return the amplitude
    if (momenta_to_integrate.empty()) {
        return amp; 
    }

    //return net energy of all momenta 
    auto get_net_energy = [](const std::array<FourVec,D>& v) {
        double E = 0.; 
        for (const auto& p : v) E += p.get(0); 
        return E; 
    }; 

    const double starting_energy = get_net_energy(momenta); 

    long long int n_accepted=0; 

    long long int i_report=0; 
    long long int step_n_accepted=0.; 
    long double step_amp_sum=0.;  

    for (long long int step=0; step<n_steps; step++) {

        amp_sum += (long double)amp; 

        double new_amp = 1.; 

        bool valid_update{false}; 
        
        auto new_momenta = momenta; 

        //go down the line, and suggest updates to each energy 
        for (int ind : momenta_to_integrate) {
            
            const auto& inp = momenta_inputs[ind]; 
            
            auto& pi = new_momenta[ind];
            pi = Nudge(
                pi, 
                Gaus()*inp.momenta_scan_amplitude, 
                Gaus()*inp.momenta_scan_amplitude, 
                Gaus()*inp.momenta_scan_amplitude
            );    
        }
        double dE = get_net_energy(new_momenta) - starting_energy;

        auto& pi = new_momenta[momenta_to_integrate.back()];
        
        //if this is the last momenta, we have to re-scale this particle's momenta to enforce energy conservation 
        double m2 = pi.norm2(); 
        double target_energy = pi[0] - dE; 
        
        //reject this update, if we can't decrease this particle's energy enough to conserve momentum 
        if (target_energy < m2) continue; 
        double p_mag = std::sqrt( target_energy*target_energy - m2 );

        double norm3 = std::sqrt( pi[0]*pi[0] - m2 ); 
        pi[0] = target_energy; 
        pi[1] *= p_mag/norm3; 
        pi[2] *= p_mag/norm3; 
        pi[3] *= p_mag/norm3; 

        //now, calculate the amplitude 
        new_amp = fcn(new_momenta); 

        if (is_nan(new_amp)) continue; 
        
        if ( Rndm() < new_amp/amp )  {
            
            //std::printf(" accepted\n");
            amp = new_amp; 
            momenta = new_momenta; 
            ++n_accepted;
            ++step_n_accepted; 
        } else {
            //std::printf(" rejected\n"); 
        }
        

        
        if (++i_report >= steps_between_reports) {
            printf(
                "step %5lli/%5lli (%5.1f%%) ----------------------------------------------------\n"
                "            total          step\n"
                "   amp_sum: %-9.4e         %-9.4e\n"
                "acceptance: %5.1f%%        %5.1f%%\n"   
                "momenta: \n", 
                step, n_steps, 100.*((double)step)/((double)n_steps),     

                (double)amp_sum/((double)step),            (double)step_amp_sum/((double)i_report),
                100.*((double)n_accepted)/((double)step),  100.*((double)step_n_accepted)/((double)i_report)
            );
            for (int i=0; i<D; i++) {
                std::cout << " " << i << " " << momenta[i] << std::endl; 
            }
            i_report=0; 
            step_amp_sum=0;
            step_n_accepted=0; 
        } else {
            step_amp_sum += amp; 
        }
    
    }
    
    return amp / ((double)n_steps);
}


#endif



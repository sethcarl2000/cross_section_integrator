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

namespace CSInteg {
    
    /// @return 'true' if arg is nan
    template<typename T> inline bool is_nan(T x) { return x!=x; };
}

struct FourVecParameters {
    bool is_fixed; //if 'true', then it will not be integrated over
    FourVec val; //starting value 
    double momenta_scan_amplitude; //amplitude of the momentum to scan
    double mass; //mass of the particle (MeV/c^2)
};

namespace Setting {
    enum Bit : int {
        kNone = 0,
        kConserveEnergy     = 1 << 0, //enforce external energy conservation for each update
        kConserveMomenta    = 1 << 1, //enforce external momentum conservation for each update 
        kAutoAdjustScan     = 1 << 2, //automatically adjust the scan rate to keep the acceptance prob around ~50% 
        kVerbose            = 1 << 3  //if true, then report output
    };
}

//integrates a square matrix element over external momenta
template<int D> long double CrossSectionIntegrator(
    std::function<long double(const std::array<FourVec,D>&)> fcn, 
    std::array<FourVecParameters,D> momenta_inputs, 
    long long int n_steps, 
    long long int steps_between_reports,
    const int mode = Setting::kConserveEnergy | Setting::kVerbose, 
    long long int n_thermalization_steps = 0 //number of steps before measurements start
)
{       
    using CSInteg::is_nan;

    const double kNaN = std::numeric_limits<double>::quiet_NaN(); 

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
    if (mode & Setting::kVerbose) std::cout << " total momenta to integrate over: " << momenta_to_integrate.size() << std::endl; 

    //try to perform monte-carlo integration of variables    
    long double amp_sum  = 0.; 
    long double amp2_sum = 0.; 

    //get the starting value of the amplitude 
    long double amp = fcn(momenta);

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

    //multiply amplitude by integral measure
    for (int i=0; i<D; i++) amp *= measure[i](momenta[i]);

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
    long double step_amp2_sum=0.; 

    long double steps_amp_average=0.; 
    long double steps_amp2_average=0.; 
    int i_reported=0; 

    for (long long int step=-n_thermalization_steps; step<n_steps; step++) {

        if (step >= 0) {
            amp_sum  += amp; 
            amp2_sum += amp*amp;   
        }

        long double new_amp = 1.; 

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
        for (int i=0; i<D; i++) new_amp *= measure[i](momenta[i]);


        if (is_nan(new_amp)) continue; 
        
        if ( Rndm() < new_amp/amp )  {
            
            //std::printf(" accepted\n");
            amp = new_amp; 
            momenta = new_momenta; 
            if (step >=0 ) ++n_accepted;
            ++step_n_accepted; 
        } else {
            //std::printf(" rejected\n"); 
        }
        

        if ((mode & Setting::kVerbose) && ++i_report >= steps_between_reports) {
            
            if (step >= 0) ++i_reported; 

            long double d_i_reported = (long double)i_reported;             
            long double d_i_report   = (long double)i_report; 
            long double d_step       = (long double)step; 

            if (step >= 0) {
                steps_amp_average  += step_amp_sum/d_i_report; 
                steps_amp2_average += step_amp2_sum/d_i_report; 
            }

            //compute the variance of the amplitude 
            long double step_rel_variance = (step_amp2_sum/d_i_report) - std::pow(step_amp_sum/d_i_report,2);
            long double variance          = (steps_amp2_average/d_i_reported) - std::pow(steps_amp_average/d_i_reported,2);

            //get the relative variance
            step_rel_variance   = std::sqrt(step_rel_variance); 
            //variance            = std::sqrt(variance);

            //step acceptance prob. 
            long double step_accept_p = ((long double)step_n_accepted)/d_i_report; 
            
            printf(
                "step %5lli/%5lli (%5.1f%%) ----------------------------------------------------\n"
                "            total                           step\n"
                "   amp_sum: %10.4Le +/- %-8.2Le             %10.4Le +/- %-8.2Le\n"
                "acceptance: %5.1Lf%%                        %5.1Lf%%\n"   
                "momenta: \n", 
                step, n_steps, 100.*((double)step)/((double)n_steps),     

                (long double)amp_sum/d_step,            std::sqrt(variance),         
                (long double)step_amp_sum/d_i_report,   step_rel_variance,
                (long double)100.*n_accepted/d_step,    (long double)100.*step_n_accepted/d_i_report
            );
            for (int i=0; i<D; i++) {
                std::cout << " " << i << " " << momenta[i] << std::endl; 
            }
            
            //number of times a report has been made
            //number of steps since the last report 
            i_report=0; 
            
            step_amp_sum=0;
            step_amp2_sum=0; 
            step_n_accepted=0; 

            //adjust the scan rate 
            if (mode & Setting::kAutoAdjustScan) {
                for (int ind : momenta_to_integrate) {
                    auto& scan_amp = momenta_inputs[ind].momenta_scan_amplitude;
                    scan_amp *= 1.  -  (0.5 - step_accept_p);
                }
            }

        } else {
        
            step_amp_sum  += amp; 
            step_amp2_sum += amp*amp; 
        
        }
    
    }
    
    return amp_sum / ((long double)n_steps);
}


#endif



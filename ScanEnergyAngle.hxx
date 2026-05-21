#ifndef ScanEnergyAngle_H
#define ScanEnergyAngle_H

//project headers
#include "FourVec.hxx"
#include "PolarFourVec.hxx"
#include "Bound.hxx"
#include "IntegSettings.hxx"
 
//root headers
#include <TH2D.h> 

//stdlib headers
#include <cmath> 
#include <stdio.h>
#include <functional> 
#include <random> 
#include <array> 
#include <iostream> 
#include <limits>

namespace SEInteg {
    
    /// @return 'true' if arg is nan
    template<typename T> inline bool is_nan(T x) { return x!=x; };

    constexpr double pi = 3.1415926536; 
}

struct FVScanParameters {
    bool is_fixed; //if 'true', then it will not be integrated over
    PolarFourVec val; //starting value 
    double mass; //mass of the particle (MeV/c^2)
    Bound<double> cos_theta_bound{ -1., +1. }; //bounds for theta (angle between incident and outgoing)
    Bound<double> energy_bound{}; //bounds for energy 
    double scan_amplitude{0.01}; //amplitude of the momentum to scan
};

/// @brief creates histogram of energy-angle distribution. 
/// @tparam D number of momenta inputs for a matrix element 
/// @param M2 input squared matrix element 
/// @param momenta_inputs momentum inputs
/// @param outputs output histograms
/// @param n_steps number of sampling steps
/// @param steps_between_reports interval between report / update steps
/// @param setting integration settings
/// @param n_thermalization_steps number of 'thermalization' steps before integration begins
/// @return measured amplitude
template<int D> long double ScanEnergyAngle(
    std::function<long double(const std::array<FourVec,D>&)> M2, 
    std::array<FVScanParameters,D> momenta_inputs, 
    std::vector<TH2D>& outputs,
    long long int n_steps, 
    long long int steps_between_reports,
    const int setting = Setting::kConserveEnergy | Setting::kVerbose, 
    long long int n_thermalization_steps = 0 //number of steps before measurements start
)
{
    double scan_rate = 0.01; 

    using SEInteg::is_nan, SEInteg::pi; 

    static const double kNaN = std::numeric_limits<double>::quiet_NaN(); 

    // a wrapper for our matrix element which takes arguments in the form of cartesian four-vectors (but we work in polar form here). 
    auto M2_polar = [&M2](const std::array<PolarFourVec,D>& P_polar) 
    {
        long double measure = 1.; 
        std::array<FourVec,D> P_FV;
        for (int i=0; i<D; i++) {
            P_FV[i] = static_cast<FourVec>{P_polar[i]}
            measure *= 0.5/P_polar[i].energy; 
        }
        return measure * M2(P_FV); 
    };
    //____________________________________________________________________________________________________
   
    //initialize random-number generators
    std::mt19937 twister; 
    std::uniform_real_distribution<double> runiform{0., 1.}; 
    auto Rndm = std::bind(runiform, twister); 

    /// @return uniform dist. random val on interval [x0, x1)
    auto Rndm_range = [&Rndm](double x0, double x1) { return x0 + (x1-x0)*Rndm(); };

    std::normal_distribution<double> rgaus{0., 1.};
    auto Gaus = std::bind(rgaus, twister); 

    //the momenta we will vary
    std::array<PolarFourVec,D> momenta; 

    //momenta to integrate (rest are fixed)
    std::vector<int> momenta_to_integrate; 

    for (int i=0; i<D; i++) {

        const auto& inp = momenta_inputs[i];

        momenta[i] = inp.val;
        if (!inp.is_fixed) { momenta_to_integrate.push_back(i); }
    }
    if (setting & Setting::kVerbose) std::cout << " total momenta to integrate over: " << momenta_to_integrate.size() << std::endl; 

    /// @return net energy of all momenta (ingoing and outgoing)
    auto get_net_energy = [](const std::array<PolarFourVec,D>& v) {
        double E{0.};
        for (const auto& pi : v) E += pi.energy; 
        return E; 
    }; 

    const double starting_energy = get_net_energy(momenta); 

    //try to perform monte-carlo integration of variables    
    long double amp_sum  = 0.; 
    long double amp2_sum = 0.; 

    //get the starting value of the amplitude 
    long double amp = M2_polar(momenta);

    if (is_nan(amp) || amp < 0.) {
        std::cerr << "<"<<__func__<<">: Fatal error: ~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "   amplitude is negative / NaN. amp=" << amp << "\n"
        "   external momenta: \n"; 
        for (int i=0; i<D; i++) {
            std::cerr << "["<<i<<"]: " << momenta[i] << "\n"; 
        }
        std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
        return kNaN; 
    }

    //if all momenta are fixed, then just return the amplitude
    if (momenta_to_integrate.empty()) {
        return amp; 
    }

    //multiply amplitude by integral measure
    for (int i=0; i<D; i++) amp *= measure[i](momenta[i]);

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

            double energy_scan_range = inp.energy_bound.span()*scan_rate; 
            
            pi.energy = inp.energy_bound.enforce( 
                pi.energy + Rndm_range(-energy_scan_range, +energy_scan_range) 
            );

            double cos_theta_scan_range = inp.cos_theta_bound.span()*scan_rate; 

            pi.cos_theta = inp.cos_theta_bound.enforce( 
                pi.cos_theta + Rndm_range(-cos_theta_scan_range, +cos_theta_scan_range)
            ); 
            
            double phi_scan_range = 2.*SEInteg::pi*scan_rate/5.; 

            pi.phi += Rndm_range( -phi_scan_range, +phi_scan_range ); 
            if (pi.phi < -SEInteg::pi) pi.phi += +2.*SEInteg::pi; 
            if (pi.phi > +SEInteg::pi) pi.phi += -2.*SEInteg::pi; 
        }
        // Since we just randomly adjusted several energies, 
        // we will shift the energy of all momenta  by this value, to enforce energy conservation.  
        double energy_adjust = (starting_energy - get_net_energy(new_momenta))/((double)D);

        for (int ind : momenta_to_integrate) {
            auto& pi = new_momenta[ind]; 
            const auto& inp = momenta_inputs[ind]; 
            pi.energy = inp.energy_bound.enforce( pi.energy + energy_adjust );  
        }
        // re-measure the energy to make sure that it matches. 
        // if it does not match, that means that, when we re-scaled all our (mutable) energies
        // to enforce energy conservation, one of our scalable energies hit a hard boundry set by 
        // the user; thus, this is an illegal update. 
        if (std::fabs(get_net_energy(new_momenta) - starting_energy) > 1.e-8) continue; 

        //now, calculate the amplitude 
        new_amp = M2_polar(new_momenta);     
        
        if (is_nan(new_amp)) continue; 
        
        if ( Rndm() < new_amp/amp )  {
            
            //std::printf(" accepted\n");
            amp = new_amp; 
            momenta = new_momenta; 
            if (step >=0 ) ++n_accepted;
            ++step_n_accepted; 

        } else {

        }
        

        if ((setting & Setting::kVerbose) && ++i_report >= steps_between_reports) {
            
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
            if (setting & Setting::kAutoAdjustScan) {
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



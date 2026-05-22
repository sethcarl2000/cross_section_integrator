#ifndef EstimateDifferentialCS_HXX
#define EstimateDifferentialCS_HXX

//project headers
#include "FourVec.hxx"
#include "PolarFourVec.hxx"
#include "Bound.hxx"
#include "IntegSettings.hxx"
#include "processes.hxx"
 
//root headers
#include <TH2D.h> 

//gsl headers
#include <gsl/gsl_monte.h> 
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

//stdlib headers
#include <cmath> 
#include <stdio.h>
#include <functional> 
#include <random> 
#include <array> 
#include <iostream> 
#include <limits>
#include <vector> 
#include <stdexcept> 
#include <sstream> 
#include <string> 
#include <algorithm>

namespace EDCS {
    constexpr double pi = 3.1415926536; 

    //set the energy scale 
    constexpr double MeV    = 1.; 
    // atomic mass unit (in MeV)
    constexpr double dalton = 931.5*MeV; 

    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN(); 

    //this struct wraps a M2_cartiesian matrix element, and accepts args as polar four-vectors
    template<int D> struct M2_polar_wrapper {

        processes::M2_fcn<D> M2_cartesian; 

        inline double operator()(const std::array<PolarFourVec,D>& v_polar) const {
            std::array<FourVec,D> v_cartesian; 
            int i=0; 
            for (const auto& p_polar : v_polar) v_cartesian[i++] = static_cast<FourVec>(p_polar);
            return (double)M2_cartesian(v_cartesian);
        }
    };
}; 

/// @brief Use GSL monte-carlo integrators to estimate the differential cross sections 
/// @tparam D number of momenta in matrix element 
/// @param M2_cartesian function representing square matrix element (taking momentum inputs in cartesian form)
/// @param P0 incoming momenta (single)
/// @param P1 single outgoing momenta to fix 
/// @param P1_ind index of outgoing momenta to fix  
/// @param n_steps number of integration steps to attempt
/// @param setting setting bit for integrator
/// @return differential cross section:  dSigma / dE * dCos (cos, E)
template<int D> double EstimateDifferentialCS(
    const processes::M2_fcn<D>& fcn_M2_cartesian, 
    const PolarFourVec& P0, 
    const PolarFourVec& P1,
    const int P1_ind,  
    const std::vector<double> spectator_mass,
    long int n_steps, 
    const double target_mass = 183.8*EDCS::dalton, 
    int setting = Setting::kVerbose
)
{
    using EDCS::pi, EDCS::kNaN; 

    static_assert(D >= 2, "D must be at least 2."); 
    
    if (2 + (int)spectator_mass.size() != D) {
        std::ostringstream oss; 
        oss << "in <EstimateDifferentialCS<"<<D<<">>: list of spectator masses provided ("<<spectator_mass.size()<<")"
        " is incorrect size. must be D-2"; 
        throw std::logic_error(oss.str()); 
        return kNaN; 
    }

    //first, we need to write a wrapper function for the M2_fcn, which accepts *polar* four-vectors. 
    // NOTE: it is assumed that P[0] (the first momentum) for this matrix element is the input momentum! 
    
    const double E_beam = P0.energy; 

    //energy of outgoing particle (to measure)
    const double E1 = P1.energy; 

    const double m2 = P1.mass2;

    if (E1 > E_beam) {
        std::ostringstream oss; 
        oss << "in <EstimateDifferentialCS<"<<D<<">>: energy of fixed outgoing particle (" << E1 << " MeV) "
        " is greater than incoming particle energy (" << E_beam << " MeV)"; 
        throw std::logic_error(oss.str()); 
        return kNaN; 
    }   

    if (P1_ind < 1) {
        std::ostringstream oss; 
        oss << "in <EstimateDifferentialCS<"<<D<<">>: index of fixed outgoing momenta = "<<P1_ind<<", is invalid; must be >0 (it cannot be the incident momenta!)"; 
        throw std::logic_error(oss.str()); 
        return kNaN; 
    }


    constexpr double twopi = 2.*pi; 

    //these are the parts of the measure that are constants. 
    // to see where this and other formulas come from (in terms of the sqaure matrix element),
    // see Peskin  Ch 4., in particular eq. 4.79. 
    double prefactor = twopi*std::sqrt( E1*E1 - m2 )/( 4.*E_beam*target_mass * std::pow(twopi,3*D) * 2. ); 

    std::array<PolarFourVec,D> momenta{ P0, P1 };

    //wrap our matrix element (accepting cartesian inputs) so that it will accept polar-form momenta inputs
    EDCS::M2_polar_wrapper<D> M2_polar{fcn_M2_cartesian};

    //in this case there are no momenta to integrate over. return the differential CS.  
    if (D==2) {
        return prefactor * M2_polar(momenta); 
    }
    
    //this is a tricky piece. How many degrees of freedom are we integrating over? 
    //first, we have to check how many 'spectator' degrees of freedom we have. 
    const int n_spectator_momenta = D-2; 
    if (setting & Setting::kVerbose) {
        printf("n. spectator momenta: %i\n", n_spectator_momenta);
    }


    //this is because the integration (in polar-coords, about the z-axis) has 3 DoF for each spectator momenta, 
    // but we must take one away, because the *last* spectator momenta has it's energy determined by whatever we need
    // to make energy conservation work. in other words, we let the energies of all spectator momenta run free, but 
    // the energy of the last spectator momenta is determined by overall energy conservation.  
    const int integral_DoF = n_spectator_momenta*3 - 1; 

    std::vector<Bound<double>> specmom_energy_bound{};
    std::vector<int> specmom_ind{}; 

    
    //bounds of the hypercube we will integrate in. They are layed out like this: 
    // [0] - P2 dCos2   [-1,+1]
    // [1] -    dPhi2   [-pi,+pi]
    // [2] -    dW      [m2, max E2] 
    // [3] - P3 dCos3   [-1,+1]
    // [4] -    dPhi3   [-pi,+pi]
    // [5] -    dW      [m3, max E3]
    // [6] - P4 dCos3   [-1,+1]
    // [7] -    dPhi3   [-pi,+pi] 
    //     <---- NO integration over energy for last 'free' momenta. 
     
    double integ_vars_lo[integral_DoF];
    double integ_vars_hi[integral_DoF];

    //this is the last spectator momentum's mass. this is the most we can reduce it's energy to! 
    const double last_spectator_mass = spectator_mass.back(); 

    if (n_spectator_momenta > 1) {
        specmom_ind.reserve(n_spectator_momenta);
        specmom_energy_bound.reserve(n_spectator_momenta);
        
        double max_energy = E_beam - E1; 
        int i_spec=0; 
        for (int i=1; i<D; i++) { if (i!=P1_ind) max_energy += -spectator_mass[i_spec++]; }
        
        i_spec=0; 
        for (int i=1; i<D; i++) {

            if (i==P1_ind) continue; 

            specmom_ind.push_back(i);
            
            integ_vars_lo[3*i_spec + 0] = -1.;
            integ_vars_hi[3*i_spec + 0] = +1.;   

            integ_vars_lo[3*i_spec + 1] = -pi;
            integ_vars_hi[3*i_spec + 1] = +pi;  

            double m = spectator_mass[i_spec]; 

            if (i_spec < n_spectator_momenta-1) {
                integ_vars_lo[3*i_spec + 2] = m;
                integ_vars_hi[3*i_spec + 2] = m + max_energy;  
#ifdef EDCS_DEBUG
                std::printf(
                    "Spectator momenta %i (momentum %i): ~~~~~~~~~~~~~~~~~~~\n"
                    "   cos(theta): [ %+5.3f, %+5.3f ]\n"
                    "   phi:        [ %+5.3f, %+5.3f ]\n"
                    "   energy:     [ %+5.3f, %+5.3f ]\n",
                    i_spec, i, 
                    integ_vars_lo[3*i_spec + 0], integ_vars_hi[3*i_spec + 0], 
                    integ_vars_lo[3*i_spec + 1], integ_vars_hi[3*i_spec + 1], 
                    integ_vars_lo[3*i_spec + 2], integ_vars_hi[3*i_spec + 2] 
                );
#endif
            } else {
#ifdef EDCS_DEBUG
                std::printf(
                    "Spectator momenta %i (momentum %i): ~~~~~~~~~~~~~~~~~~~\n"
                    "   cos(theta): [ %+5.3f, %+5.3f ]\n"
                    "   phi:        [ %+5.3f, %+5.3f ]\n",
                    i_spec, i, 
                    integ_vars_lo[3*i_spec + 0], integ_vars_hi[3*i_spec + 0], 
                    integ_vars_lo[3*i_spec + 1], integ_vars_hi[3*i_spec + 1]
                );
#endif
            }
            i_spec++; 
        }
    }

#ifdef EDCS_DEBUG
    std::printf("total momenta: %i, spectator momenta: %i\n");
#endif

    //now, we're ready to assemble the function that will be passed to the gsl function. 
    auto my_MC_integrand = [&specmom_energy_bound, &specmom_ind, &momenta, &M2_polar, P1_ind, P0, P1, n_spectator_momenta, E_beam, prefactor]
        (double *X, size_t dim, void *params)
    {   
        std::array<PolarFourVec,D> P; P[0] = P0; 

        double E_balance = P0.energy - P1.energy; 
        
        //we're going to 'unfold' the integrator's inputs into our momenta  
        int i_spec = 0; 
        for (int i=1; i<D; i++) {

            if (i==P1_ind) {
                P[i] = P1;  
            } else {
                auto& pi = P[i]; 
                pi.cos_theta    = X[3*i_spec + 0]; 
                pi.phi          = X[3*i_spec + 1];
                if (i_spec < n_spectator_momenta-1) {
                    pi.energy   = X[3*i_spec + 2];
                    E_balance  -= pi.energy; 
                } 
                ++i_spec;
            }
        }

        //now, lets pick the energy of the final spec. momentum 
        auto& p_pivot = P[specmom_ind.back()]; 
        if (E_balance < std::sqrt( p_pivot.mass2 )) return 0.; 
        
        p_pivot.energy = E_balance; 

#ifdef EDCS_DEBUG
        std::printf(
            "matrix element will be evaluated.. ~~~~~~~~~~~~~~~~~~~~\n"  
        );
        std::cout << "specmom indicies:"; 
        for (int ind : specmom_ind) std::cout << " " << ind; 
        std::cout << std::endl; 

        double E = 0.; 
        for (int i=0; i<D; i++) {
            
            bool is_spectator = 
                !(std::find(specmom_ind.begin(), specmom_ind.end(), i) == specmom_ind.end());
            
            const auto& p = P[i]; 

            E += p.energy; 

            std::printf(
                "Momentum P[%1i] (%s) ----------------------\n"
                "   cos(theta): %- 5.3f \n"
                "   phi:        %- 5.3f \n"
                "   energy:     %- 5.3f \n"
                "   mass^2      %- 5.3f \n",
                i, (is_spectator ? "spectator" : "fixed"),
                p.cos_theta, 
                p.phi, 
                p.energy, 
                p.mass2
            );
        }
        std::printf(
            "energy balance: %.3f%%\n", 
            100.*(E - P0.energy)/P0.energy
        );
#endif
        return prefactor * M2_polar(P); 
    };

    

    double x[integral_DoF] = {
        0.99, pi/2., 500., 
        0.98,-pi/2. 
    };

    return my_MC_integrand(x, integral_DoF, nullptr); 
} 

#endif
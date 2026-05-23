
#include "DifferentialCSIntegrator.hxx"

//project headers
#include "FourVec.hxx"
#include "PolarFourVec.hxx"
#include "Bound.hxx"
#include "IntegSettings.hxx"
#include "processes.hxx"
 
//root headers
#include <Math/AdaptiveIntegratorMultiDim.h>
#include <Math/AllIntegrationTypes.h>
#include <Math/IntegratorMultiDim.h>
#include <Math/IntegratorOptions.h>
#include <Math/Functor.h>

//gsl headers
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h> 
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>

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
#include <random> 

namespace {
    double monte_f(double* X, size_t dim, void* fcn_wrapper_void)
    {
        const std::function<double(double*)>* fcn = (std::function<double(double*)>*)(fcn_wrapper_void);       
        return (*fcn)(X);
    }
}

IntegrationResult DifferentialCSIntegrator::Integrate(
    const PolarFourVec& P0, 
    const PolarFourVec& P1,
    int P1_ind,
    const std::vector<double> spectator_mass
) 
{
    const double kNaN = some_numbers::NaN<double>; 

    using 
        some_numbers::pi, 
        some_numbers::twopi;

    IntegrationResult result; 
    
    const int D = fExpr.get_n_inputs() - 2; 

    if (2 + (int)spectator_mass.size() != D) {
        std::ostringstream oss; 
        oss << "in <DifferentialCDIntegrator::"<<__func__<<">: list of spectator masses provided ("<<spectator_mass.size()<<")"
        " is incorrect size. must be D-2"; 
        throw std::logic_error(oss.str()); 
        return (result | ResultStatus::kError); 
    }

    //first, we need to write a wrapper function for the M2_fcn, which accepts *polar* four-vectors. 
    // NOTE: it is assumed that P[0] (the first momentum) for this matrix element is the input momentum! 
    
    const double E_beam = P0.energy; 

    //energy of outgoing particle (to measure)
    const double E1 = P1.energy; 

    const double m2 = P1.mass2;

    if (E1 > E_beam) {
        std::ostringstream oss; 
        oss << "in <DifferentialCDIntegrator::"<<__func__<<">: energy of fixed outgoing particle (" << E1 << " MeV) "
        " is greater than incoming particle energy (" << E_beam << " MeV)"; 
        throw std::logic_error(oss.str()); 
        return (result | ResultStatus::kError); 
    }   

    if (P1_ind < 1) {
        std::ostringstream oss; 
        oss << "in <DifferentialCDIntegrator::"<<__func__<<">: index of fixed outgoing momenta = "<<P1_ind<<", is invalid; must be >0 (it cannot be the incident momenta!)"; 
        throw std::logic_error(oss.str()); 
        return (result | ResultStatus::kError); 
    }

    const double target_mass = some_numbers::mZ; 

    //these are the parts of the measure that are constants. 
    // to see where this and other formulas come from (in terms of the sqaure matrix element),
    // see Peskin  Ch 4., in particular eq. 4.79. 
    double prefactor = twopi*std::sqrt( E1*E1 - m2 )/( 4.*E_beam*target_mass * std::pow(twopi,3*D) * 2. ); 

    
    std::vector<PolarFourVec> momenta; momenta.reserve(D);
    momenta.push_back(P0);
        
    //in this case there are no momenta to integrate over. return the differential CS.  
    if (D==2) {
        momenta.push_back(P1);
        return IntegrationResult{ 
            .val = prefactor * fExpr(momenta),
            .error = 0., 
            ResultStatus::kSuccess 
        }; 
    }
    
    //this is a tricky piece. How many degrees of freedom are we integrating over? 
    //first, we have to check how many 'spectator' degrees of freedom we have. 
    const int n_spectator_momenta = D-2; 
    if (fOptions & Setting::kVerbose) {
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

    if (n_spectator_momenta > 0) {
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
            } 
#ifdef EDCS_DEBUG
            else {
                std::printf(
                    "Spectator momenta %i (momentum %i): ~~~~~~~~~~~~~~~~~~~\n"
                    "   cos(theta): [ %+5.3f, %+5.3f ]\n"
                    "   phi:        [ %+5.3f, %+5.3f ]\n",
                    i_spec, i, 
                    integ_vars_lo[3*i_spec + 0], integ_vars_hi[3*i_spec + 0], 
                    integ_vars_lo[3*i_spec + 1], integ_vars_hi[3*i_spec + 1]
                );
            }
#endif
            
            i_spec++; 
        }
    }

#ifdef EDCS_DEBUG
    std::printf("total momenta: %i, spectator momenta: %i\n", D, n_spectator_momenta);
#endif

    //now, we're ready to assemble the function that will be passed to the gsl function. 
    auto my_MC_integrand = [&specmom_energy_bound, &specmom_ind, &momenta, this, D, P1_ind, P0, P1, n_spectator_momenta, E_beam, prefactor]
        (double *X)
    {   
        std::vector<PolarFourVec> P; P.reserve(D);
        P[0] = P0; 

        double E_balance = P0.energy - P1.energy; 
        
        //we're going to 'unfold' the integrator's inputs into our momenta  
        int i_spec = 0; 
        for (int i=1; i<D; i++) {

            if (i==P1_ind) {
                P.push_back(P1);   
            } else {
                PolarFourVec pi; 
                pi.cos_theta    = X[3*i_spec + 0]; 
                pi.phi          = X[3*i_spec + 1];
                if (i_spec < n_spectator_momenta-1) {
                    pi.energy   = X[3*i_spec + 2];
                    E_balance  -= pi.energy; 
                } 
                P.push_back(pi);
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
        double M2 = fExpr(P);
        if (M2 < 0. || M2 != M2) return 0.; 
        
        return prefactor * M2; 
    };

    std::function<double(double*)> my_MC_integrand_f{my_MC_integrand};

    //same function as above, but accepts const ptr of arguments
    std::function<double(const double*)> my_MC_integrand_fconst = [&my_MC_integrand_f, integral_DoF](const double* X)
    {
        double xx[integral_DoF];
        std::copy( X, X+integral_DoF, xx );
        return my_MC_integrand_f(xx);
    };

    gsl_monte_function gsl_M2_func;

    gsl_M2_func.f      = monte_f;    
    gsl_M2_func.dim    = integral_DoF; 
    gsl_M2_func.params = &my_MC_integrand_f; 
    
#ifdef EDCS_DEBUG 
    std::cout << "setting up gsl rng environment...\n";
#endif
    const gsl_rng_type *type;
    gsl_rng *rng;

    gsl_rng_env_setup();
     
    type = gsl_rng_default;
    rng = gsl_rng_alloc(type);

    std::random_device rd; 
    gsl_rng_set(rng, rd());

    double result{kNaN}, error{kNaN}; 

    switch (fStrategy) {
        
        //____________________________________________________________________________________________________________
        case Setting::kVEGAS : {
            //initialize the state of the monte-carlo integrator
            gsl_monte_vegas_state *vegas_state = gsl_monte_vegas_alloc(integral_DoF); 

            gsl_monte_vegas_params vegas_params; 
            gsl_monte_vegas_params_get(vegas_state, &vegas_params); 

            vegas_params.alpha = 2.; 
            vegas_params.iterations = 10; 
            vegas_params.ostream = stdout; 
            vegas_params.verbose = -1;     

            gsl_monte_vegas_params_set(vegas_state, &vegas_params); 

            //now, we're ready to do the integration. 
            gsl_monte_vegas_integrate(
                &gsl_M2_func, 
                integ_vars_lo, integ_vars_hi, integral_DoF, 
                (size_t)fMaxCalls, 
                rng, vegas_state, 
                &result.val, &result.error
            );

            result = result | ResultStatus::kSuccess; 

            //free memory 
            gsl_monte_vegas_free(vegas_state);
            break; 
        }  

        //____________________________________________________________________________________________________________
        case Setting::kMISER : {
            //initialize the state of the monte-carlo integrator
            gsl_monte_miser_state *miser_state = gsl_monte_miser_alloc(integral_DoF); 

            gsl_monte_miser_params miser_params; 
            gsl_monte_miser_params_get(miser_state, &miser_params); 

            miser_params.alpha = 2.; 
            miser_params.min_calls = 16.*integral_DoF;
            miser_params.min_calls_per_bisection = 16 * miser_params.min_calls; 
            
            gsl_monte_miser_params_set(miser_state, &miser_params); 

            //now, we're ready to do the integration. 
            gsl_monte_miser_integrate(
                &gsl_M2_func, 
                integ_vars_lo, integ_vars_hi, integral_DoF, 
                (size_t)fMaxCalls, 
                rng, miser_state, 
                &result.val, &result.error 
            );

            result = result | ResultStatus::kSuccess; 

            //free memory 
            gsl_monte_miser_free(miser_state);
            break; 
        }

        //____________________________________________________________________________________________________________
        case Setting::kPLAIN : {
            //initialize the state of the monte-carlo integrator
            gsl_monte_plain_state *plain_state = gsl_monte_plain_alloc(integral_DoF); 

            //now, we're ready to do the integration. 
            gsl_monte_plain_integrate(
                &gsl_M2_func, 
                integ_vars_lo, integ_vars_hi, integral_DoF, 
                (size_t)fMaxCalls, 
                rng, plain_state, 
                &result.val, &result.error 
            );

            result = result | ResultStatus::kSuccess; 

            //free memory 
            gsl_monte_plain_free(plain_state);
            break; 
        }

        //____________________________________________________________________________________________________________
        case Setting::kADAPTIVE : {

            ROOT::Math::IntegratorMultiDim integ(ROOT::Math::IntegrationMultiDim::kADAPTIVE);

            auto options = integ.Options(); 
            options.SetNCalls((unsigned int)fMaxCalls);
            integ.SetOptions(options);

            integ.SetRelTolerance(fRelTolerance);

            result.val = integ.Integral(my_MC_integrand_fconst, integral_DoF, integ_vars_lo, integ_vars_hi);
            result.error = integ.Error(); 
            
            //check the integrator status
            int status = integ.Status(); 

            switch (status) {
                case 0  : result.flag = ResultStatus::kSuccess; break; 
                case 1  : result.flag = ResultStatus::kToleranceNotAchieved | ResultStatus::kError; break; 
                default : result.flag = ResultStatus::kError; 
            }

            if ((fOptions & Setting::kVerbose) && (result.flag & ResultStatus::kToleranceNotAchieved)) {
                std::printf("Desired rel. tolerance not achieved within max_calls = %li\n", fMaxCalls); 
            }
            
            break; 
        }
        
        default : {
            std::cerr << "unsupported integration techniuqe.\n" << std::endl; return (result | ResultStatus::kError); 
        }
    }
    

#ifdef EDCS_DEBUG 
    std::cout << "monte-carlo integration starting...\n";
#endif

    gsl_rng_free(rng);

    if (fOptions & Setting::kVerbose)
        printf("done with integration. result: %.5e +/- %.5e\n", result, error);

    return result; 
}
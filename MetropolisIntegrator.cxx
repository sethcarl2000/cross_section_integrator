#include "MetropolisIntegrator.hxx"
#include "Bound.hxx"

#include <cstdio> 

//__________________________________________________________________________________________
MetropolisIntegrator::MetropolisIntegrator(const std::function<double(const double*)>& _fcn, int _DoF)
    : fFcn{_fcn}, fDoF{_DoF}
{
    //initialize random number generator
    std::random_device rd; 
    fTwister = std::mt19937{rd()};   

    fState = std::vector<double>(fDoF, 0.);
}
//__________________________________________________________________________________________
IntegrationResult MetropolisIntegrator::Integrate(const double* x_min, const double *x_max)
{
    //find the biggest value with a greedy search 
    std::vector<Bound<double>> bounds; bounds.reserve(fDoF);
    for (int i=0; i<fDoF; i++) {

        bounds.emplace_back(x_min[i], x_max[i]);
        fState[i] = (x_min[i] + x_max[i])/2.;
    }

    double scan_rate =1.;

    auto new_state = fState; 

    fAmplitude = Amplitude(fState); 

    for (unsigned int i_greedy=0; i_greedy<fNSteps_greedy; i_greedy++) {
    
        for (int i=0; i<fDoF; i++) {

            auto& x = fState[i];

            const auto& bound = bounds[i];
            new_state[i] = bound.enforce( x + bound.span()*scan_rate*RandRange(-0.5, +0.5) ); 
        }

        double new_amplitude;
        if ((new_amplitude = Amplitude(new_state)) > fAmplitude) {
            fState = new_state; 
            fAmplitude = new_amplitude; 
            std::printf("accepted. ~~~~~~~~~~~~~~~~");
        }
        std::printf("rejected.");
        std::printf(" new amplitude: %e\n", fAmplitude);
        
        scan_rate *= std::pow(0.001, 1./((double)fNSteps_greedy));
    }

    return IntegrationResult{};
}
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
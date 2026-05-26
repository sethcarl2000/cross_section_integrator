    #ifndef MetropolisIntegrator_HXX
#define MetropolisIntegrator_HXX

#include <cmath> 
#include <vector> 
#include <functional> 
#include <random> 
#include "IntegrationResult.hxx"

class MetropolisIntegrator {
private: 

    std::uniform_real_distribution<double> fRand{0., 1.};
    std::mt19937 fTwister; 

    inline double Rand() { return fRand(fTwister); }
    inline double RandRange(double min, double max) { return min + (max-min)*Rand(); }

    const std::function<double(const double*)>& fFcn; 

    inline double Amplitude(const std::vector<double>& X) { return fFcn(X.data()); }

    int fDoF;
    std::vector<double> fState; 
    double fAmplitude; 

    double fScanAmplitude{0.1}; 

    unsigned int fNSteps_greedy{(unsigned int)100e3};
    unsigned int fNSteps_integration{(unsigned int)100e6};  

    unsigned int fNUpdates{100}; 

    double fTarget_probability{0.20}; 

public: 

    MetropolisIntegrator(const std::function<double(const double*)>& _fcn, int _DoF); 

    // number of integration steps
    void SetIntegrationSteps(unsigned int steps) { fNSteps_integration=steps; }
    
    // number of 'greedy' steps (initial search for maximum, beofore integration)
    void SetGreedySteps(unsigned int steps) { fNSteps_greedy=steps; }

    // goal for fraction of updates accepted (higher is less aggresive walk, lower is more aggresive walk)
    void SetTargetProbability(double target_prob) { fTarget_probability=target_prob; }

    // number of times to update the search amplitude 
    void SetNUpdates(unsigned int n_updates) { fNUpdates=n_updates; }

    // set the starting scan-amplitude
    void SetScanAmplitude(double scan_amplitude) { fScanAmplitude=scan_amplitude; }

    IntegrationResult Integrate(const double *X_min, const double *X_max);
   
};


#endif
#ifndef DifferentialCSIntegrator_HXX
#define DifferentialCSIntegrator_HXX

#include "processes.hxx"
#include "IntegSettings.hxx"
#include "IntegrationResult.hxx"

class DifferentialCSIntegrator {

    processes::M2_expr fExpr; 

    long int fMaxCalls{3e7};    

    int fOptions{Setting::kVerbose}; 

    double fRelTolerance; 

    Setting::MCStrategy fStrategy{Setting::kADAPTIVE}; 

private: 

    DifferentialCSIntegrator(const processes::M2_expr& _expr) : fExpr{_expr} {}; 

    void SetProcess(const processes::M2_expr& _expr) { fExpr = _expr; }

    void SetMaxCalls(long int n_steps) { fMaxCalls = n_steps; }

    void SetIntegrationStrategy(Setting::MCStrategy strat) { fStrategy=strat; }

    void SetOptions(int options) { fOptions=options; }

    void SetRelTolerance(double tolerance) { fRelTolerance=tolerance; }

    IntegrationResult Integrate(
        const PolarFourVec& P0, 
        const PolarFourVec& P1,
        int P1_ind,
        const std::vector<double> spectator_mass
    );

};

#endif
#include "processes.hxx"
#include "some_numbers.hxx"

namespace processes 
{
    
M2_expr Factory::elastic()
{   
    std::function<double(const std::vector<FourVec>&)> _fcn = [](const std::vector<FourVec>& P)
    { 
        using namespace some_numbers; 
    
        static const double prefactor = Square(4*pi*e0*Z)/4.;

        const auto& P0 = P[0]; 
        const auto& P1 = P[1]; 
        
        auto Q = P1 - P0; 

        //electron / positron mass
        static const double m = me, m2 = m*m;

        const double  Q2 = Q*Q; 

        const double  eP0 = P0.get(0); 
        
        double M2 = 8.*eP0*eP0 - 2.*Q2;

        return prefactor * M2 / (Q2*Q2); 
    };

    return M2_expr(2, _fcn);
} 
    
}
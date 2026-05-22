#include "processes.hxx"

namespace processes 
{
    
    long double bh_photoproduction(const std::array<FourVec,3>& P)
    {   
        static const long double prefactor = 1817.91495662; // (4pi*Z*e^2)^2 / 4

        const auto& P0 = P[0]; 
        const auto& P1 = P[1]; 
        const auto& K  = P[2]; 

        auto Q = (P1 + K) - P0; 

        //electron / positron mass
        static const long double m = 0.501, m2 = m*m, m4 = m2*m2, m6 = m4*m2, m8 = m6*m2; 

        const long double  P0P1 = P0*P1;
        const long double  P0Q  = P0*Q; 
        const long double  P1Q  = P1*Q;  
        const long double  Q2   = Q*Q; 

        const long double  eP0 = P0.get(0); 
        const long double  eP1 = P1.get(0);

        long double M2 = (-16*P0Q*P1Q)/Square(2*P1Q - Q2) - (16*eP0*eP1*Q2)/Square(2*P1Q - Q2) - (16*m2*Q2)/Square(2*P1Q - Q2) + (8*P0P1*Q2)/Square(2*P1Q - Q2) - (16*P0Q*P1Q)/Square(2*P0Q + Q2) - (256*Square(eP1)*m2*Square(P0Q))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (128*Square(eP1)*P0P1*Square(P0Q))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (512*eP0*eP1*m2*P0Q*P1Q)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) - (256*eP0*eP1*P0P1*P0Q*P1Q)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) - (256*Square(eP0)*m2*Square(P1Q))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (128*Square(eP0)*P0P1*Square(P1Q))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (64*Square(eP1)*Square(P0Q))/((2*P1Q - Q2)*Square(2*P0Q + Q2)) - (128*eP0*eP1*P0Q*P1Q)/((2*P1Q - Q2)*Square(2*P0Q + Q2)) + (64*Square(eP0)*Square(P1Q))/((2*P1Q - Q2)*Square(2*P0Q + Q2)) - (16*eP0*eP1*Q2)/Square(2*P0Q + Q2) - (16*m2*Q2)/Square(2*P0Q + Q2) + (8*P0P1*Q2)/Square(2*P0Q + Q2) - (256*eP0*eP1*m2*P0Q*Q2)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) - (256*Square(eP1)*m2*P0Q*Q2)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (128*eP0*eP1*P0P1*P0Q*Q2)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (128*Square(eP1)*P0P1*P0Q*Q2)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (256*Square(eP0)*m2*P1Q*Q2)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (256*eP0*eP1*m2*P1Q*Q2)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) - (128*Square(eP0)*P0P1*P1Q*Q2)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) - (128*eP0*eP1*P0P1*P1Q*Q2)/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (32*eP0*eP1*P0Q*Q2)/((2*P1Q - Q2)*Square(2*P0Q + Q2)) + (32*Square(eP1)*P0Q*Q2)/((2*P1Q - Q2)*Square(2*P0Q + Q2)) - (32*Square(eP0)*P1Q*Q2)/((2*P1Q - Q2)*Square(2*P0Q + Q2)) - (32*eP0*eP1*P1Q*Q2)/((2*P1Q - Q2)*Square(2*P0Q + Q2)) - (64*Square(eP0)*m2*Square(Q2))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) - (128*eP0*eP1*m2*Square(Q2))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) - (64*Square(eP1)*m2*Square(Q2))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (32*Square(eP0)*P0P1*Square(Q2))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (64*eP0*eP1*P0P1*Square(Q2))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) + (32*Square(eP1)*P0P1*Square(Q2))/(Square(2*P1Q - Q2)*Square(2*P0Q + Q2)) - (64*Square(eP1)*Square(P0Q))/(Square(2*P1Q - Q2)*(2*P0Q + Q2)) + (128*eP0*eP1*P0Q*P1Q)/(Square(2*P1Q - Q2)*(2*P0Q + Q2)) - (64*Square(eP0)*Square(P1Q))/(Square(2*P1Q - Q2)*(2*P0Q + Q2)) - (32*eP0*eP1*P0Q*Q2)/(Square(2*P1Q - Q2)*(2*P0Q + Q2)) - (32*Square(eP1)*P0Q*Q2)/(Square(2*P1Q - Q2)*(2*P0Q + Q2)) + (32*Square(eP0)*P1Q*Q2)/(Square(2*P1Q - Q2)*(2*P0Q + Q2)) + (32*eP0*eP1*P1Q*Q2)/(Square(2*P1Q - Q2)*(2*P0Q + Q2)) + (16*P0P1*Q2)/((2*P1Q - Q2)*(2*P0Q + Q2));
    
        return prefactor * M2 / (Q2*Q2); 
    } 
    
}
#ifndef processes_HXX
#define processes_HXX

#include "FourVec.hxx"
#include "CrossSectionIntegrator.hxx"
#include <array>
#include <cmath> 
#include <functional> 

namespace processes 
{
    //invalid for i<2!; 
    double Power(double x, int i) { do { x *= x; } while (--i > 1); return x; }
    
    /// @brief trident process e-(P0)   -->   e-(P1)  e-(Pm)   e+(Pp)
    /// @param P[0] P0 incoming beam electron 
    /// @param P[1] P1 outgoing beam electron 
    /// @param P[2] Pm generated pair electron (outgoing)
    /// @param P[3] Pp generated pair positron (outgoing) 
    /// @return Squre amplitude for given momentum 
    double trident(const std::array<FourVec,4>& P);
};

#endif
#ifndef processes_HXX
#define processes_HXX

#include "FourVec.hxx"
#include <array>
#include <cmath> 
#include <functional> 

namespace processes 
{
    template<int D> using M2_fcn = std::function<long double(const std::array<FourVec,D>&)>; 

    /// @brief trident process e-(P0)   -->   e-(P1)  e-(Pm)   e+(Pp)
    /// @param P[0] P0 incoming beam electron 
    /// @param P[1] P1 outgoing beam electron 
    /// @param P[2] Pm generated pair electron (outgoing)
    /// @param P[3] Pp generated pair positron (outgoing) 
    /// @return Squre amplitude for given momentum 
    long double trident(const std::array<FourVec,4>& P);
};

#endif
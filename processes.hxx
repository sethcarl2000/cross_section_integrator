#ifndef processes_HXX
#define processes_HXX

#include "FourVec.hxx"
#include <array>
#include <cmath> 
#include <functional> 

namespace processes 
{
    inline long double Square(long double _x) { return _x*_x; }

    template<int D> using M2_fcn = std::function<long double(const std::array<FourVec,D>&)>; 

    /// @brief trident process e-(P0)   -->   e-(P1)  e-(Pm)   e+(Pp)
    /// @param P[0] P0 incoming beam electron 
    /// @param P[1] P1 outgoing beam electron 
    /// @param P[2] Pm generated pair electron (outgoing)
    /// @param P[3] Pp generated pair positron (outgoing) 
    /// @return Squre amplitude for given momentum (including averging over external spins)
    long double trident(const std::array<FourVec,4>& P);

    /// @brief photoproduction from bethe-heitler process: e-(P0)   -->     e-(P1)  gamma(K)
    /// @param P[0] P0 incoming electron momentum
    /// @param P[1] P1 outgoing electron momentum 
    /// @param P[2] K  outgoing photon momentum 
    /// @return Square amplitude (including averging over spin-sums)
    long double bh_photoproduction(const std::array<FourVec,3>& P);
};

#endif
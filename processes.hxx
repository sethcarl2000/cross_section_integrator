#ifndef processes_HXX
#define processes_HXX

#include "FourVec.hxx"
#include "PolarFourVec.hxx"
#include <array>
#include <vector> 
#include <cmath> 
#include <functional> 

namespace processes 
{
    inline long double Square(long double _x) { return _x*_x; }

    class M2_expr {
    private: 
        size_t n_inputs =1.; 
        std::function<double(const std::vector<FourVec>&)> fM2_cartiesian;
    public: 

        M2_expr(size_t n_inp, const std::function<double(const std::vector<FourVec>&)>& fM2_cart); 

        int get_n_inputs() const { return n_inputs; }

        //operator on cartesian coordinates
        double operator()(const std::vector<FourVec>&) const; 

        //operator on polar coordinates
        double operator()(const std::vector<PolarFourVec>&) const; 
    };


    /// @brief trident process e-(P0)   -->   e-(P1)  e-(Pm)   e+(Pp)
    /// @param P[0] P0 incoming beam electron 
    /// @param P[1] P1 outgoing beam electron 
    /// @param P[2] Pm generated pair electron (outgoing)
    /// @param P[3] Pp generated pair positron (outgoing) 
    /// @return Squre amplitude for given momentum (including averging over external spins)
    M2_expr trident();

    /// @brief photoproduction from bethe-heitler process: e-(P0)   -->     e-(P1)  gamma(K)
    /// @param P[0] P0 incoming electron momentum
    /// @param P[1] P1 outgoing electron momentum 
    /// @param P[2] K  outgoing photon momentum 
    /// @return Square amplitude (including averging over spin-sums)
    M2_expr bh_photoproduction();
};

#endif
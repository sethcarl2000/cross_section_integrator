#ifndef some_numbers_HXX
#define some_numbers_HXX

#include <limits> 

//some numbers
namespace some_numbers 
{
    constexpr double pi = 3.1415926535, twopi = 2.*pi; 

    //value of the electric charge, in natural units
    constexpr double e0 = 0.302822120768; 

    //charge of tungsten nucleus 
    constexpr double Z = 74.; 

    //some units
    constexpr double MeV = 1.; 
    constexpr double dalton = 931.5*MeV; 

    // electron mass 
    constexpr double me = 0.511*MeV, me2 = me*me; 

    // tungsten nucleus mass 
    constexpr double mZ = 183.8*dalton; 

    // nan
    template <typename T> T NaN = std::numeric_limits<T>::quiet_NaN();  
};   

#endif
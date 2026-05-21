#ifndef PolarFourVec_HXX
#define PolarFourVec_HXX

#include "FourVec.hxx"
#include <cmath> 
#include <string>
#include <iostream>  

//an alternative represetntation of FourVec; 
struct PolarFourVec {

    double cos_theta, phi, energy, mass2; 

    PolarFourVec() {}; 

    //conversion operator from FourVec to PolarFourVec 
    explicit inline PolarFourVec(const FourVec& v) 
    {
        //cosine of theta
        cos_theta = v.z() / std::sqrt( v.x()*v.x() + v.y()*v.y() + v.z()*v.z() ) ; 

        //phi (azimuthal angle)
        if (v.x()==v.y()==0.) { //handle the case of an exactly polar vector
            phi = 0.;
        } else {    
            phi = std::acos( v.x() / std::sqrt( v.x()*v.x() + v.y()*v.y() ) ) * std::signbit(v.y()) ? -1 : +1; 
        }
        
        energy  =  v.t();
        mass2   =  v.norm2();
    }; 

    //conversion from PolarFourVec to FourVec
    
    //explicit conversion from FourVec to PolarFourVec
    explicit inline operator FourVec() const {
        double p = std::sqrt( energy*energy - mass2 );
        double sin_theta  = std::sqrt( 1. - cos_theta*cos_theta ); 
        return FourVec{{
            energy, 
            p*std::cos(phi)*sin_theta, 
            p*std::sin(phi)*sin_theta, 
            p*cos_theta
        }};
    }; 
};
 
/// @return string formatted with PolarFourVec info 
std::string String(const PolarFourVec&); 

std::ostream& operator<<(std::ostream& os, const PolarFourVec& v);

//check to make sure the polar four-vec is trivially copy- & move-constructable
static_assert(std::is_trivially_copy_constructible_v<PolarFourVec>, "PolarFourVec is not trivially copy-constructable"); 
static_assert(std::is_trivially_move_constructible_v<PolarFourVec>, "PolarFourVec is not trivially move-constructable"); 

#endif
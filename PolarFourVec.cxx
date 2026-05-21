#include "PolarFourVec.hxx"
#include "FourVec.hxx"

#include <string> 

//_________________________________________________________________________________
std::string String(const PolarFourVec& v) 
{
    char buff[50]; 
    std::sprintf(buff, "(p0 = % 5.4f, cos(theta) = % 5.4f, phi = % 5.4f rad, m^2 = % 5.4f)", v.energy, v.cos_theta, v.phi, v.mass2);
    return std::string{buff}; 
}
//_________________________________________________________________________________
std::ostream& operator<<(std::ostream& os, const PolarFourVec& v)
{
    return os << String(v); 
}
//_________________________________________________________________________________
//_________________________________________________________________________________
//_________________________________________________________________________________
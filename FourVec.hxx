#ifndef FourVec_HXX
#define FourVec_HXX

#include <array>
#include <cmath> 
#include <string> 
#include <iostream> 

struct FourVec {
    std::array<double,4> P;

    inline double& operator[](int i) { return P[i]; }
    inline double get(int i) const { return P[i]; }

    inline double operator*(const FourVec& rhs) const {
        return this->get(0)*rhs.get(0) - (this->get(1)*rhs.get(1) + this->get(2)*rhs.get(2) + this->get(3)*rhs.get(3)); 
    }

    //norm2
    inline double norm2() const { return (*this)*(*this); }
    inline double norm()  const { return std::sqrt(norm2()); }
};  
 
FourVec Nudge(FourVec v, double x, double y, double z) {
    double m2 = v.norm2(); 
    v[1] += x; v[2] += y; v[3] += z; 
    v[0] = std::sqrt( m2 + v[1]*v[1] + v[2]*v[2] + v[3]*v[3] );
    return v; 
} 

inline FourVec operator+(const FourVec& lhs, const FourVec& rhs) { 
    return {
        lhs.get(0) + rhs.get(0), 
        lhs.get(1) + rhs.get(1), 
        lhs.get(2) + rhs.get(2), 
        lhs.get(3) + rhs.get(3)
    }; 
}

inline FourVec operator-(const FourVec& lhs, const FourVec& rhs) { 
    return {
        lhs.get(0) - rhs.get(0), 
        lhs.get(1) - rhs.get(1), 
        lhs.get(2) - rhs.get(2), 
        lhs.get(3) - rhs.get(3)
    }; 
}


std::string String(const FourVec& v) {
    char buff[50]; 
    std::sprintf(buff, "(% 5.1f, % 5.1f % 5.1f % 5.1f)", v.P[0], v.P[1],v.P[2],v.P[3]);
    return std::string{buff}; 
}

std::ostream& operator<<(std::ostream& os, const FourVec& v) {
    return os << String(v); 
}

static_assert(std::is_trivially_copy_constructible_v<FourVec>, "FourVec is not trivially copy-constructable"); 
static_assert(std::is_trivially_move_constructible_v<FourVec>, "FourVec is not trivially move-constructable"); 


#endif
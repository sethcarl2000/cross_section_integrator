#ifndef Bound_HXX
#define Bound_HXX

template<typename T> struct Bound { 
    T min, max; 
    inline bool operator()(const T& x) const { return x > min && x < max; }
    inline T span() const { return max-min; }

    //if a given value is outside the bound, re-enforce it to the nearest bound
    inline T enforce(T val) const {
        if (val < min) return min;
        if (val > max) return max; 
        return val;  
    }
};

#endif
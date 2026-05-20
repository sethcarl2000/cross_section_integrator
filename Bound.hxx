#ifndef Bound_HXX
#define Bound_HXX

template<typename T> struct Bound { 
    T min, max; 
    inline bool operator()(const T& x) const { return x > min && x < max; }
};

#endif
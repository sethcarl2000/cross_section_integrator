#ifndef IntegrationResult_HXX
#define IntegrationResult_HXX

#include "some_numbers.hxx"

namespace ResultStatus {
    enum Flag : int {
        kNull    = 0 << 0,
        kSuccess = 1 << 0, 
        kError   = 1 << 1,
        // when using the 'adaptive' method, this bit is used to represent when the max calls is reached before 
        // the desired tolerance is achieved. 
        kToleranceNotAchieved = 1 << 2 
    }; 
};

struct IntegrationResult {

    double val{some_numbers::NaN<double>}, error{some_numbers::NaN<double>}; 
    int flag{ResultStatus::kNull}; 

    //implicit conversion to double
    operator double() const { return val; }

    //explicit conversion to bool 
    explicit operator bool() const { return (flag & ResultStatus::kSuccess); }

    inline IntegrationResult operator|(int rhs) const {
        return IntegrationResult{
            .val = val, 
            .error = error, 
            .flag = flag | rhs
        };  
    }
};

#endif
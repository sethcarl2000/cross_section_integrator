#include "processes.hxx"
#include "some_numbers.hxx"
#include <stdexcept>
#include <sstream> 
#include <string> 

namespace {
    const double kNaN = some_numbers::NaN<double>; 

    std::function<double(std::vector<FourVec>&)> null_fcn = [](const std::vector<FourVec>&){ return kNaN; }; 
};

namespace processes 
{

//____________________________________________________________________________________________
M2_expr::M2_expr(size_t n_inp, const std::function<double(const std::vector<FourVec>&)>& M2_cart)
    : n_inputs{n_inp}, fM2_cartesian{M2_cart}
{
    //default constructor 
};
//____________________________________________________________________________________________
double M2_expr::operator()(const std::vector<FourVec>& P) const
{
    if (P.size() == n_inputs) {

        return fM2_cartesian(P);
    }

    std::ostringstream oss;
    oss << "in <M2_expr::"<<__func__<<">: invalid number of inupts. expected " << n_inputs << ", got " << P.size(); 
    throw std::logic_error(oss.str());

    return kNaN; 
}
//____________________________________________________________________________________________
double M2_expr::operator()(const std::vector<PolarFourVec>& P) const
{
    if (P.size() == n_inputs) {

        //we have to convert these polar four-vecs to cartiesan four-vecs
        std::vector<FourVec> P_cartesian; P_cartesian.reserve(n_inputs); 
        for (const auto& P_polar : P) {

            //use the explicit conversion operator
            P_cartesian.emplace_back((FourVec)P_polar); 
        }
        return fM2_cartesian(P_cartesian);
    }

    std::ostringstream oss;
    oss << "in <M2_expr::"<<__func__<<">: invalid number of inupts. expected " << n_inputs << ", got " << P.size(); 
    throw std::logic_error(oss.str());
    
    return kNaN; 
}
//____________________________________________________________________________________________
//____________________________________________________________________________________________
//____________________________________________________________________________________________
//____________________________________________________________________________________________
};
/******************************************************************************
 * @file Tensor.cpp
 * @brief Stream output for the Tensor class
 *****************************************************************************/

#include "Tensor.hpp"

#include <iomanip>
#include <ostream>


std::ostream& operator<<(std::ostream& os, const Tensor& T)
{
    // save the current format
    auto flags = os.flags();
    auto prec = os.precision();

    // change the format for tensor output
    os << std::fixed << std::setprecision(6);

    // tensor output: row-major bracketed rows
    os  << '(' << T.xx() << ", " << T.xy() << ", " << T.xz() << "; "
        << T.yx() << ", " << T.yy() << ", " << T.yz() << "; "
        << T.zx() << ", " << T.zy() << ", " << T.zz() << ')';

    // restore the original format
    os.flags(flags);
    os.precision(prec);

    return os;
}

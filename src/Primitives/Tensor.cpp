/******************************************************************************
 * @file Tensor.cpp
 * @brief Stream output for the Tensor class
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Tensor.h"

// Standard library headers
#include <iomanip>
#include <ostream>

// *************************** Non-Member Functions ***************************

std::ostream& operator<<(std::ostream& os, const Tensor& T)
{
    // save the current format
    const auto flags = os.flags();
    const auto prec = os.precision();

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

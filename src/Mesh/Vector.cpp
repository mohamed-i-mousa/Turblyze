/******************************************************************************
 * @file Vector.cpp
 * @brief Stream output for the Vector class
 *****************************************************************************/

#include <ostream>
#include <iomanip>

#include "Vector.hpp"


std::ostream& operator<<(std::ostream& os, const Vector& p)
{
    std::ios_base::fmtflags flags = os.flags();
    auto prec = os.precision();

    os  << std::fixed << std::setprecision(6);

    os  << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";

    os.flags(flags);
    os.precision(prec);

    return os;
}

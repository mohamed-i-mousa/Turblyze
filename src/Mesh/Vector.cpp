/******************************************************************************
 * @file Vector.cpp
 * @brief Stream output for the Vector class
 *****************************************************************************/

#include "Vector.hpp"

#include <ostream>
#include <iomanip>


std::ostream& operator<<(std::ostream& os, const Vector& p)
{
    auto flags = os.flags();
    auto prec = os.precision();

    os  << std::fixed << std::setprecision(6);

    os  << '(' << p.x() << ", " << p.y() << ", " << p.z() << ')';

    os.flags(flags);
    os.precision(prec);

    return os;
}

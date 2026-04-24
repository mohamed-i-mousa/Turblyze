/******************************************************************************
 * @file Vector.cpp
 * @brief Stream output for the Vector class
 *****************************************************************************/

#include "Vector.h"

#include <ostream>
#include <iomanip>


std::ostream& operator<<(std::ostream& os, const Vector& p)
{
    // save the current format 
    auto flags = os.flags();
    auto prec = os.precision();

    // change the format for vector output
    os  << std::fixed << std::setprecision(6);

    // vector output
    os  << '(' << p.x() << ", " << p.y() << ", " << p.z() << ')';

    // restore the original format
    os.flags(flags);
    os.precision(prec);

    return os;
}

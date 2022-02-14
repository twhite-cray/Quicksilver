#ifndef PHYSICAL_CONSTANTS_HH
#define PHYSICAL_CONSTANTS_HH

#include "DeclareMacro.hh"
namespace PhysicalConstants
{

static constexpr double _neutronRestMassEnergy = 9.395656981095e+2; /* MeV */
static constexpr double _pi = 3.1415926535897932;
static constexpr double _speedOfLight  = 2.99792458e+10;                // cm / s

// Constants used in math for computer science, roundoff, and other reasons
static constexpr double _tinyDouble           = 1.0e-13;
static constexpr double _smallDouble          = 1.0e-10;
static constexpr double _hugeDouble           = 1.0e+75;
//
}


#endif

#ifndef MCT_DISTANCE_INCLUDE
#define MCT_DISTANCE_INCLUDE

#include "DeclareMacro.hh"

class MC_Distance_To_Facet
{
public:
    double distance;
    int facet;
    int subfacet;
    __host__ __device__ MC_Distance_To_Facet(): distance(0.0), facet(0), subfacet(0) {}
private:
    MC_Distance_To_Facet( const MC_Distance_To_Facet& );                    // disable copy constructor
    MC_Distance_To_Facet& operator=( const MC_Distance_To_Facet& tmp );     // disable assignment operator

};

#endif

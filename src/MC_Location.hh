#ifndef MC_LOCATION_INCLUDE
#define MC_LOCATION_INCLUDE


// ToDo:  How much chaos would be caused by removing the default constructor?

#include <string>
#include "DeclareMacro.hh"

class  MC_Domain;
class  MC_Cell_State;
class  MonteCarlo;

class MC_Location
{
 public:
   int domain;
   int cell;
   int facet;

   __host__ __device__ MC_Location(int adomain, int acell, int afacet)
   : domain(adomain),
     cell(acell),
     facet(afacet)
   {}

   MC_Location()
   : domain(-1),
     cell(-1),
     facet(-1)
   {}

   const MC_Domain& get_domain(MonteCarlo *mcco) const;
};

inline bool operator==(const MC_Location& a, const MC_Location b)
{
   return
      a.domain == b.domain &&
      a.cell == b.cell &&
      a.facet == b.facet;
}


#endif

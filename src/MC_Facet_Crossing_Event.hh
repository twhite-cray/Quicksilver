#ifndef MC_FACET_CROSSING_EVENT_HH
#define MC_FACET_CROSSING_EVENT_HH

#include "MC_Tally_Event.hh"
#include "DeclareMacro.hh"

//----------------------------------------------------------------------------------------------------------------------
//  Determines whether the particle has been tracked to a facet such that it:
//    (i) enters into an adjacent cell
//   (ii) escapes across the system boundary (Vacuum BC), or
//  (iii) reflects off of the system boundary (Reflection BC).
//
//----------------------------------------------------------------------------------------------------------------------


static inline MC_Tally_Event::Enum MC_Facet_Crossing_Event(MC_Particle &mc_particle, Device &device, const int particle_index, const int maxCount, int *__restrict__ const sendCounts, MessageParticle *__restrict__ const sendParts)
{

    MC_Location location = mc_particle.Get_Location();
    const DeviceFacet &facet = device.domains[location.domain].cells[location.cell].facets[location.facet];

    if ( facet.event == MC_Subfacet_Adjacency_Event::Transit_On_Processor )
    {
        // The particle will enter into an adjacent cell.
        mc_particle.domain     = facet.adjacentDomain;
        mc_particle.cell       = facet.adjacentCell;
        mc_particle.facet      = facet.adjacentFacet;
        mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Transit_Exit;
    }
    else if ( facet.event == MC_Subfacet_Adjacency_Event::Boundary_Escape )
    {
        // The particle will escape across the system boundary.
        mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Escape;
    }
    else if ( facet.event == MC_Subfacet_Adjacency_Event::Boundary_Reflection )
    {
        // The particle will reflect off of the system boundary.
        mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Reflection;
    }
    else if ( facet.event == MC_Subfacet_Adjacency_Event::Transit_Off_Processor )
    {
        // The particle will enter into an adjacent cell on a spatial neighbor.
        // The neighboring domain is on another processor. Set domain local domain on neighbor proc
        
        mc_particle.domain     = facet.adjacentDomain;
        mc_particle.cell       = facet.adjacentCell;
        mc_particle.facet      = facet.adjacentFacet;
        mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Communication;

        device.processing[particle_index] = mc_particle;

        const int neighbor_rank = device.domains[facet.currentDomain].neighbors[facet.neighbor];
        const int offset = atomicAdd(sendCounts+neighbor_rank,1);
        sendParts[maxCount*neighbor_rank+offset] = mc_particle;
    }

    return mc_particle.last_event;
}

#endif


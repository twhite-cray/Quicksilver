#include "CycleTracking.hh"
#include "MonteCarlo.hh"
#include "MC_Base_Particle.hh"
#include "MC_Segment_Outcome.hh"
#include "CollisionEvent.hh"
#include "MC_Facet_Crossing_Event.hh"
#include "MCT.hh"
#include "DeclareMacro.hh"
#include "AtomicMacro.hh"
#include "macros.hh"
#include "qs_assert.hh"

void CycleTrackingGuts( const int numParticles, Device &device, const int maxCount, int *const sendCounts, MessageParticle *const sendParts)
{
    MC_Particle mc_particle;
    int previous = -1;
    int particle_index = 0;
    while (particle_index < numParticles) {

        if (previous != particle_index) {
          previous = particle_index;
          device.processing[particle_index].set(mc_particle);
          if (mc_particle.time_to_census <= 0) mc_particle.time_to_census += device.timeStep;
          mc_particle.energy_group = device.getEnergyGroup(mc_particle.kinetic_energy);
        }

        // Determine the outcome of a particle at the end of this segment such as:
        //
        //   (0) Undergo a collision within the current cell,
        //   (1) Cross a facet of the current cell,
        //   (2) Reach the end of the time step and enter census,
        //
        MC_Segment_Outcome_type::Enum segment_outcome = MC_Segment_Outcome(device, mc_particle);

        ATOMIC_UPDATE( device.tallies[Device::SEGMENTS] );

        mc_particle.num_segments += 1.;  /* Track the number of segments this particle has
                                            undergone this cycle on all processes. */
        switch (segment_outcome) {
        case MC_Segment_Outcome_type::Collision:
            {
              if (CollisionEvent(device, mc_particle) != MC_Collision_Event_Return::Continue_Tracking) {
                particle_index++;
              }
            }
            break;
    
        case MC_Segment_Outcome_type::Facet_Crossing:
            {
                // The particle has reached a cell facet.
                MC_Tally_Event::Enum facet_crossing_type = MC_Facet_Crossing_Event(mc_particle, device, particle_index, maxCount, sendCounts, sendParts);

                if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Transit_Exit)
                {}
                else if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Escape)
                {
                    ATOMIC_UPDATE( device.tallies[Device::ESCAPE] );
                    mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Escape;
                    particle_index++;
                }
                else if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Reflection)
                {
                    MCT_Reflect_Particle(device, mc_particle);
                }
                else
                {
                    // Enters an adjacent cell in an off-processor domain.
                    particle_index++;
                }
            }
            break;
    
        case MC_Segment_Outcome_type::Census:
            {
                const int iProcessed = __atomic_fetch_add(device.particleSizes+Device::PROCESSED,1,__ATOMIC_RELAXED);
                device.processed[iProcessed] = mc_particle;
                ATOMIC_UPDATE( device.tallies[Device::CENSUS] );
                particle_index++;
                break;
            }
            
        default:
           qs_assert(false);
           break;  // should this be an error
        }
    }
}


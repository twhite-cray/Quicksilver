#include "CycleTracking.hh"
#include "MonteCarlo.hh"
#include "ParticleVaultContainer.hh"
#include "ParticleVault.hh"
#include "MC_Segment_Outcome.hh"
#include "CollisionEvent.hh"
#include "MC_Facet_Crossing_Event.hh"
#include "MCT.hh"
#include "DeclareMacro.hh"
#include "AtomicMacro.hh"
#include "macros.hh"
#include "qs_assert.hh"

void CycleTrackingGuts( MonteCarlo *monteCarlo, int numParticles, ParticleVault *processingVault, ParticleVault *processedVault , Device &device)
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
            // The particle undergoes a collision event producing:
            //   (0) Other-than-one same-species secondary particle, or
            //   (1) Exactly one same-species secondary particle.
              if (CollisionEvent(monteCarlo->_particleVaultContainer, device, mc_particle) != MC_Collision_Event_Return::Continue_Tracking) {
                processingVault->invalidateParticle( particle_index );
                device.processing[particle_index++].species = -1;
              }
            }
            break;
    
        case MC_Segment_Outcome_type::Facet_Crossing:
            {
                // The particle has reached a cell facet.
                MC_Tally_Event::Enum facet_crossing_type = MC_Facet_Crossing_Event(mc_particle, monteCarlo, particle_index, processingVault);

                if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Transit_Exit)
                {}
                else if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Escape)
                {
                    ATOMIC_UPDATE( monteCarlo->_tallies->_balanceTask[0]._escape);
                    mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Escape;
                    mc_particle.species = -1;
                    processingVault->invalidateParticle( particle_index );
                    device.processing[particle_index++].species = -1;
                }
                else if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Reflection)
                {
                    MCT_Reflect_Particle(monteCarlo, mc_particle);
                }
                else
                {
                    // Enters an adjacent cell in an off-processor domain.
                    processingVault->invalidateParticle( particle_index );
                    device.processing[particle_index++].species = -1;
                }
            }
            break;
    
        case MC_Segment_Outcome_type::Census:
            {
                // The particle has reached the end of the time step.
                processedVault->pushParticle(mc_particle);
                const int iProcessed = __atomic_fetch_add(device.particleSizes+Device::PROCESSED,1,__ATOMIC_RELAXED);
                device.processed[iProcessed] = mc_particle;
                ATOMIC_UPDATE( monteCarlo->_tallies->_balanceTask[0]._census);
                processingVault->invalidateParticle( particle_index );
                device.processing[particle_index++].species = -1;
                break;
            }
            
        default:
           qs_assert(false);
           break;  // should this be an error
        }
    }
}


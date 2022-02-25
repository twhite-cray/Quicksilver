#ifndef COLLISION_EVENT_HH
#define COLLISION_EVENT_HH

#include "DeclareMacro.hh"
#include "Device.hh"
#include "MacroscopicCrossSection.hh"
#include "MC_Particle.hh"
#include "MC_RNG_State.hh"
#include "PhysicalConstants.hh"

//----------------------------------------------------------------------------------------------------------------------
//  Routine MC_Collision_Event determines the isotope, reaction and secondary (projectile)
//  particle characteristics for a collision event.
//
//  Return true if the particle will continue.
//----------------------------------------------------------------------------------------------------------------------

__host__ __device__ static inline void updateTrajectory( const double energy, const double angle, MC_Particle& particle )
{
    particle.kinetic_energy = energy;
    const double cosTheta = angle;
    double randomNumber = rngSample(&particle.random_number_seed);
    const double phi = 2 * 3.14159265 * randomNumber;
    const double sinPhi = sin(phi);
    const double cosPhi = cos(phi);
    const double sinTheta = sqrt((1.0 - (cosTheta*cosTheta)));
    particle.direction_cosine.Rotate3DVector(sinTheta, cosTheta, sinPhi, cosPhi);
    const double speed = (PhysicalConstants::_speedOfLight *
            sqrt((1.0 - ((PhysicalConstants::_neutronRestMassEnergy *
            PhysicalConstants::_neutronRestMassEnergy) /
            ((energy + PhysicalConstants::_neutronRestMassEnergy) *
            (energy + PhysicalConstants::_neutronRestMassEnergy))))));
    particle.velocity.x = speed * particle.direction_cosine.alpha;
    particle.velocity.y = speed * particle.direction_cosine.beta;
    particle.velocity.z = speed * particle.direction_cosine.gamma;
    randomNumber = rngSample(&particle.random_number_seed);
    particle.num_mean_free_paths = -1.0*log(randomNumber);
}


__host__ __device__ static inline bool CollisionEvent(Device &device, MC_Particle &mc_particle, __shared__ unsigned long *__restrict__ const tallies)
{
   const int globalMatIndex = device.domains[mc_particle.domain].cells[mc_particle.cell].material;

   //------------------------------------------------------------------------------------------------------------------
   //    Pick the isotope and reaction.
   //------------------------------------------------------------------------------------------------------------------
   double randomNumber = rngSample(&mc_particle.random_number_seed);
   double totalCrossSection = mc_particle.totalCrossSection;
   double currentCrossSection = totalCrossSection * randomNumber;
   int selectedIso = -1;
   int selectedUniqueNumber = -1;
   int selectedReact = -1;
   const int numIsos = device.mats[globalMatIndex].isoSize;
   
   for (int isoIndex = 0; isoIndex < numIsos && currentCrossSection >= 0; isoIndex++)
   {
      const int uniqueNumber = device.mats[globalMatIndex].isos[isoIndex].gid;
      const int numReacts = device.reactionSize;
      for (int reactIndex = 0; reactIndex < numReacts; reactIndex++)
      {
         currentCrossSection -= macroscopicCrossSection(device, reactIndex, mc_particle.domain, mc_particle.cell,
                   isoIndex, mc_particle.energy_group);
         if (currentCrossSection < 0)
         {
            selectedIso = isoIndex;
            selectedUniqueNumber = uniqueNumber;
            selectedReact = reactIndex;
            break;
         }
      }
   }
   qs_assert(selectedIso != -1);

   //------------------------------------------------------------------------------------------------------------------
   //    Do the collision.
   //------------------------------------------------------------------------------------------------------------------
   atomicFetchAdd(tallies+Device::COLLISION,1UL);
   const NuclearDataReaction::Enum reactionType = device.isotopes[selectedUniqueNumber].reactions[selectedReact+1].type;
   switch(reactionType) {
     case NuclearDataReaction::Scatter:
       { 
         atomicFetchAdd(tallies+Device::SCATTER,1UL);
         const double mat_mass = device.mats[globalMatIndex].mass;
         const double energyOut = mc_particle.kinetic_energy*(1.0-(rngSample(&mc_particle.random_number_seed)*(1.0/mat_mass)));
         const double angleOut = rngSample(&mc_particle.random_number_seed)*2.0-1.0;
         updateTrajectory( energyOut, angleOut, mc_particle);
         mc_particle.energy_group = device.getEnergyGroup(mc_particle.kinetic_energy);
         return true;
       }
       break;
     case NuclearDataReaction::Absorption:
       atomicFetchAdd(tallies+Device::ABSORB,1UL);
       return false;
       break;
     case NuclearDataReaction::Fission:
       { 
         const int nOut = int(device.nuBar+rngSample(&mc_particle.random_number_seed));
         if( nOut == 0 ) return false;
         const int extraOffset = (nOut > 1) ? atomicFetchAdd(device.particleSizes+Device::EXTRAS,nOut) : 0;
         const double ran = rngSample(&mc_particle.random_number_seed)/2.0+0.5;
         const double energyOut = 20.0*ran*ran;
         const double angleOut = rngSample(&mc_particle.random_number_seed)*2.0-1.0;
         for (int i = 1; i < nOut; i++) {
           const double ran = rngSample(&mc_particle.random_number_seed)/2.0+0.5;
           const double energyOut = 20.0*ran*ran;
           const double angleOut = rngSample(&mc_particle.random_number_seed)*2.0-1.0;
           MC_Particle secondaryParticle = mc_particle;
           secondaryParticle.random_number_seed = rngSpawn_Random_Number_Seed(&mc_particle.random_number_seed);
           secondaryParticle.identifier = secondaryParticle.random_number_seed;
           updateTrajectory( energyOut, angleOut, secondaryParticle );
           device.extras[extraOffset+i] = secondaryParticle;
         }
         updateTrajectory( energyOut, angleOut, mc_particle);

         // If a fission reaction produces secondary particles we also add the original
         // particle to the "extras" that we will handle later.  This avoids the 
         // possibility of a particle doing multiple fission reactions in a single
         // kernel invocation and overflowing the extra storage with secondary particles.
         if ( nOut > 1 ) device.extras[extraOffset] = mc_particle;

         atomicFetchAdd(tallies+Device::FISSION,1UL);
         atomicFetchAdd<unsigned long>(tallies+Device::PRODUCE,nOut);
         mc_particle.energy_group = device.getEnergyGroup(mc_particle.kinetic_energy);
         return nOut == 1;
       }
       break;
     case NuclearDataReaction::Undefined:
       abort();
   }
   abort();
   return false;
}

#endif


#include "DeclareMacro.hh"

// Forward Declaration
class ParticleVault;
class MonteCarlo;
class MC_Particle;

void CycleTrackingGuts( MonteCarlo *monteCarlo, int particle_index, ParticleVault *processingVault, ParticleVault *processedVault );

void CycleTrackingFunction( MonteCarlo *monteCarlo, MC_Particle &mc_particle, int particle_index, ParticleVault* processingVault, ParticleVault* processedVault);

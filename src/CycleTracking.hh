#include "DeclareMacro.hh"

// Forward Declaration
struct Device;
class ParticleVault;
class MonteCarlo;
class MC_Particle;

void CycleTrackingGuts( MonteCarlo *monteCarlo, int particle_index, ParticleVault *processingVault, ParticleVault *processedVault, Device &device );

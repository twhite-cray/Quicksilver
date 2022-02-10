#include "DeclareMacro.hh"

// Forward Declaration
struct Device;
class MonteCarlo;
class MC_Particle;

void CycleTrackingGuts( MonteCarlo *monteCarlo, int particle_index, Device &device );

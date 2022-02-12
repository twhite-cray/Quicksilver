#include "DeclareMacro.hh"

// Forward Declaration
struct Device;
struct MessageParticle;

void CycleTrackingGuts( int particle_index, Device &device, int maxCount, int *sendCounts, MessageParticle *sendParts );

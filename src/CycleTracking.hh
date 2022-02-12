#include "DeclareMacro.hh"

// Forward Declaration
struct Device;
struct MessageParticle;

void CycleTrackingGuts( int ipLo, int ipHi, Device &device, int maxCount, int *sendCounts, MessageParticle *sendParts );

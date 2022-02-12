#ifndef MC_FACET_CROSSING_EVENT_HH
#define MC_FACET_CROSSING_EVENT_HH

#include "MC_Tally_Event.hh"
#include "DeclareMacro.hh"

struct Device;
class MC_Particle;
struct MessageParticle;

MC_Tally_Event::Enum MC_Facet_Crossing_Event(MC_Particle &mc_particle, Device &device, int particle_index, int maxCount, int *sendCounts, MessageParticle *sendParts);

#endif


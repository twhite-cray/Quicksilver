#ifndef COLLISION_EVENT_HH
#define COLLISION_EVENT_HH

class MonteCarlo;
class MC_Particle;

#include "DeclareMacro.hh"
bool CollisionEvent(MonteCarlo* monteCarlo, MC_Particle &mc_particle, unsigned int tally_index );


#endif


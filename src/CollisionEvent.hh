#ifndef COLLISION_EVENT_HH
#define COLLISION_EVENT_HH

struct Device;
class MC_Particle;
class ParticleVaultContainer;

#include "DeclareMacro.hh"
bool CollisionEvent(ParticleVaultContainer *particleVaultContainer, Device &device, MC_Particle &mc_particle);


#endif


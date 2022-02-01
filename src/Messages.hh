#pragma once
#include "cudaUtils.hh"
#include <mpi.h>

class MC_Base_Particle;
class MonteCarlo;

struct MessageParticle {
  MessageParticle(const MC_Base_Particle &that);
  MessageParticle(const MessageParticle &that) = default;
  void set(MC_Base_Particle &that);
  double3 coordinate;
  double3 velocity;
  double kineticEnergy;
  double weight;
  double timeToCensus;
  double age;
  double numMeanFreePaths;
  double numSegments;
  uint64_t randomNumberSeed;
  uint64_t identifier;
  int numCollisions;
  int breed;
  int species;
  int domain;
  int cell;
};

struct Messages {
  Messages();
  ~Messages();
  void init(MonteCarlo &mc);
  int nMessages;
  int maxCount;
  int *counts;
  int *ranks;
  MessageParticle *recvParts;
  MPI_Request *recvReqs;
  MessageParticle *sendParts;
  MPI_Request *sendReqs;
};


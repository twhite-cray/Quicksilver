#pragma once
#include "cudaUtils.hh"
#include <mpi.h>

class MC_Base_Particle;
class MonteCarlo;

struct MessageParticle {
  MessageParticle &operator=(const MC_Base_Particle &that);
  bool operator==(const MC_Base_Particle &that) const;
  void set(MC_Base_Particle &that) const;
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
  Messages(MonteCarlo &mc);
  ~Messages();
  void init();
  void completeRecvs();
  void completeSends();
  void startRecvs();
  void startSends();
  void addParticle(const MC_Base_Particle &part, int buffer); 

  MonteCarlo &mc;
  int tag;
  MPI_Datatype mpiParticle;
  int nMessages;
  int maxCount;
  int *counts;
  int *ranks;
  MessageParticle *recvParts;
  MPI_Request *recvReqs;
  MPI_Status *recvStats;
  MessageParticle *sendParts;
  MPI_Request *sendReqs;
};


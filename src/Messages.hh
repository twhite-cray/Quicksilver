#pragma once
#include "cudaUtils.hh"
#include "MC_Base_Particle.hh"
#include <mpi.h>

class MonteCarlo;

struct MessageParticle {

  MessageParticle &operator=(const MC_Base_Particle &that)
  {
    coordinate.x = that.coordinate.x;
    coordinate.y = that.coordinate.y;
    coordinate.z = that.coordinate.z;
    velocity.x = that.velocity.x;
    velocity.y = that.velocity.y;
    velocity.z = that.velocity.z;
    kineticEnergy = that.kinetic_energy;
    weight = that.weight;
    timeToCensus = that.time_to_census;
    age = that.age;
    numMeanFreePaths = that.num_mean_free_paths;
    numSegments = that.num_segments;
    randomNumberSeed = that.random_number_seed;
    identifier = that.identifier;
    numCollisions = that.num_collisions;
    breed = that.breed;
    domain = that.domain;
    cell = that.cell;
    return *this;
  }

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
  int *ranks;
  int *recvCounts;
  MessageParticle *recvParts;
  MPI_Request *recvReqs;
  MPI_Status *recvStats;
  int *sendCounts;
  MessageParticle *sendParts;
  MPI_Request *sendReqs;
};


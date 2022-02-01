#include "Messages.hh"

#include "cudaUtils.hh"
#include "MC_Base_Particle.hh"
#include "MC_Particle_Buffer.hh"

MessageParticle::MessageParticle(const MC_Base_Particle &that):
  coordinate{that.coordinate.x,that.coordinate.y,that.coordinate.z},
  velocity{that.velocity.x,that.velocity.y,that.velocity.z},
  kineticEnergy(that.kinetic_energy),
  weight(that.weight),
  timeToCensus(that.time_to_census),
  age(that.age),
  numMeanFreePaths(that.num_mean_free_paths),
  numSegments(that.num_segments),
  randomNumberSeed(that.random_number_seed),
  identifier(that.identifier),
  numCollisions(that.num_collisions),
  breed(that.breed),
  species(that.species),
  domain(that.domain),
  cell(that.cell)
{}

void MessageParticle::set(MC_Base_Particle &that)
{
  that.coordinate.x = coordinate.x;
  that.coordinate.y = coordinate.y;
  that.coordinate.z = coordinate.z;
  that.velocity.x = velocity.x;
  that.velocity.y = velocity.y;
  that.velocity.z = velocity.z;
  that.kinetic_energy = kineticEnergy;
  that.weight = weight;
  that.time_to_census = timeToCensus;
  that.age = age;
  that.num_mean_free_paths = numMeanFreePaths;
  that.num_segments = numSegments;
  that.random_number_seed = randomNumberSeed;
  that.identifier = identifier;
  that.last_event = MC_Tally_Event::Facet_Crossing_Communication;
  that.num_collisions = numCollisions;
  that.breed = breed;
  that.species = species;
  that.domain = domain;
  that.cell = cell;
}

Messages::Messages():
  nMessages(0),
  maxCount(0),
  counts(nullptr),
  ranks(nullptr),
  recvParts(nullptr),
  recvReqs(nullptr),
  sendParts(nullptr),
  sendReqs(nullptr)
{}

Messages::~Messages()
{
  nMessages = 0;
  maxCount = 0;
  CHECK(hipHostFree(counts)); counts = nullptr;
  delete [] ranks;
  CHECK(hipHostFree(recvParts)); recvParts = nullptr;
  delete [] recvReqs;
  CHECK(hipHostFree(sendParts)); sendParts = nullptr;
  delete [] sendReqs;
}

void Messages::init(MonteCarlo &mc)
{
  nMessages = mc.particle_buffer->num_buffers;
  maxCount = mc.particle_buffer->buffer_size;
  CHECK(hipHostMalloc(&counts,nMessages*sizeof(*counts)));
  ranks = new int[nMessages];
  recvReqs = new MPI_Request[nMessages];
  sendReqs = new MPI_Request[nMessages];
  for (int i = 0; i < nMessages; i++) {
    ranks[i] = MPI_PROC_NULL;
    recvReqs[i] = sendReqs[i] = MPI_REQUEST_NULL;
  }
  assert(nMessages == mc.particle_buffer->processor_buffer_map.size());
  for (const auto &pair : mc.particle_buffer->processor_buffer_map) {
    assert(pair.second < nMessages);
    ranks[pair.second] = pair.first;
  }
  const size_t msgBytes = sizeof(MessageParticle)*nMessages*maxCount;
  CHECK(hipHostMalloc(&recvParts,msgBytes));
  CHECK(hipHostMalloc(&sendParts,msgBytes));
}


#include "Messages.hh"

#include "cudaUtils.hh"
#include "MC_Base_Particle.hh"
#include "MC_Particle_Buffer.hh"

static int uniqueTag = 33;

MessageParticle &MessageParticle::operator=(const MC_Base_Particle &that)
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
  species = that.species;
  domain = that.domain;
  cell = that.cell;
  return *this;
}

bool MessageParticle::operator==(const MC_Base_Particle &that) const
{
  return ((coordinate.x == that.coordinate.x) &&
      (coordinate.y == that.coordinate.y) &&
      (coordinate.z == that.coordinate.z) &&
      (velocity.x == that.velocity.x) &&
      (velocity.y == that.velocity.y) &&
      (velocity.z == that.velocity.z) &&
      (kineticEnergy == that.kinetic_energy) &&
      (weight == that.weight) &&
      (timeToCensus == that.time_to_census) &&
      (age == that.age) &&
      (numMeanFreePaths == that.num_mean_free_paths) &&
      (numSegments == that.num_segments) &&
      (randomNumberSeed == that.random_number_seed) &&
      (identifier == that.identifier) &&
      (that.last_event == MC_Tally_Event::Facet_Crossing_Communication) &&
      (numCollisions == that.num_collisions) &&
      (breed == that.breed) &&
      (species == that.species) &&
      (domain == that.domain) &&
      (cell == that.cell));
}

void MessageParticle::set(MC_Base_Particle &that) const
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

Messages::Messages(MonteCarlo &mc_):
  mc(mc_),
  tag(uniqueTag++),
  mpiParticle(MPI_DATATYPE_NULL),
  nMessages(0),
  maxCount(0),
  counts(nullptr),
  ranks(nullptr),
  recvParts(nullptr),
  recvReqs(nullptr),
  recvStats(nullptr),
  sendParts(nullptr),
  sendReqs(nullptr)
{
  MPI_Type_contiguous(sizeof(MessageParticle),MPI_BYTE,&mpiParticle);
  MPI_Type_commit(&mpiParticle);
}

Messages::~Messages()
{
  delete [] sendReqs;
  CHECK(hipHostFree(sendParts)); sendParts = nullptr;
  delete [] recvStats;
  delete [] recvReqs;
  CHECK(hipHostFree(recvParts)); recvParts = nullptr;
  delete [] ranks;
  CHECK(hipHostFree(counts)); counts = nullptr;
  maxCount = 0;
  nMessages = 0;
  MPI_Type_free(&mpiParticle);
}

void Messages::init()
{
  if (nMessages) {
    assert(maxCount == mc.particle_buffer->buffer_size);
    assert(nMessages == mc.particle_buffer->processor_buffer_map.size());
    for (const auto &pair : mc.particle_buffer->processor_buffer_map) {
      assert(pair.second < nMessages);
      assert(ranks[pair.second] == pair.first);
    }
    return;
  }

  nMessages = mc.particle_buffer->num_buffers;
  maxCount = mc.particle_buffer->buffer_size;
  CHECK(hipHostMalloc(&counts,nMessages*sizeof(*counts)));
  ranks = new int[nMessages];
  recvReqs = new MPI_Request[nMessages];
  recvStats = new MPI_Status[nMessages];
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

void Messages::addParticle(const MC_Base_Particle &part, const int buffer)
{
  assert(buffer < nMessages);
  const int i = counts[buffer]++;
  assert(i < maxCount);
  sendParts[buffer*maxCount+i] = part;
}

void Messages::completeRecvs()
{
  MPI_Waitall(nMessages,recvReqs,recvStats);
}

void Messages::completeSends()
{
  MPI_Waitall(nMessages,sendReqs,MPI_STATUSES_IGNORE);
  for (int i = 0; i < nMessages; i++) counts[i] = 0;
}

void Messages::startRecvs()
{
  for (int i =0, offset = 0; i < nMessages; i++, offset += maxCount) {
    assert(recvReqs[i] == MPI_REQUEST_NULL);
    MPI_Irecv(recvParts+offset,maxCount,mpiParticle,ranks[i],tag,MPI_COMM_WORLD,recvReqs+i);
  }
}

void Messages::startSends()
{
  for (int i =0, offset = 0; i < nMessages; i++, offset += maxCount) {
    assert(sendReqs[i] == MPI_REQUEST_NULL);
    assert(counts[i] <= maxCount);
    MPI_Isend(sendParts+offset,counts[i],mpiParticle,ranks[i],tag,MPI_COMM_WORLD,sendReqs+i);
  }
}

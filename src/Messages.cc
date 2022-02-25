#include "Messages.hh"

#include "cudaUtils.hh"
#include "MC_Base_Particle.hh"
#include "MC_Particle_Buffer.hh"
#include "MonteCarlo.hh"

static int uniqueTag = 33;

bool MessageParticle::operator==(const MC_Base_Particle &that) const
{
  assert(coordinate.x == that.coordinate.x);
  assert(coordinate.y == that.coordinate.y);
  assert(coordinate.z == that.coordinate.z);
  assert(velocity.x == that.velocity.x);
  assert(velocity.y == that.velocity.y);
  assert(velocity.z == that.velocity.z);
  assert(kineticEnergy == that.kinetic_energy);
  assert(weight == that.weight);
  assert(timeToCensus == that.time_to_census);
  assert(age == that.age);
  assert(numMeanFreePaths == that.num_mean_free_paths);
  assert(numSegments == that.num_segments);
  assert(randomNumberSeed == that.random_number_seed);
  assert(identifier == that.identifier);
  assert(that.last_event == MC_Tally_Event::Facet_Crossing_Communication);
  assert(numCollisions == that.num_collisions);
  assert(breed == that.breed);
  assert(domain == that.domain);
  assert(cell == that.cell);
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
  that.domain = domain;
  that.cell = cell;
}

Messages::Messages():
  tag(uniqueTag++),
  mpiParticle(MPI_DATATYPE_NULL),
  nMessages(0),
  maxCount(0),
  ranks(nullptr),
  recvCounts(nullptr),
  recvParts(nullptr),
  recvReqs(nullptr),
  recvStats(nullptr),
  sendCounts(nullptr),
  sendParts(nullptr),
  sendReqs(nullptr)
{
  MPI_Type_contiguous(sizeof(MessageParticle),MPI_BYTE,&mpiParticle);
  MPI_Type_commit(&mpiParticle);
}

Messages::~Messages()
{
  CHECK(hipEventDestroy(event));
  delete [] sendReqs;
  CHECK(hipFree(sendParts)); sendParts = nullptr;
  CHECK(hipFree(sendCounts)); sendCounts = nullptr;
  delete [] recvStats;
  delete [] recvReqs;
  CHECK(hipFree(recvParts)); recvParts = nullptr;
  CHECK(hipFree(recvCounts)); recvCounts = nullptr;
  delete [] ranks;
  maxCount = 0;
  nMessages = 0;
  MPI_Type_free(&mpiParticle);
}

void Messages::init(MonteCarlo &mc)
{
  static constexpr int BUFFER_FRACTION = 8;
  if (nMessages) {
    assert(nMessages == mc.particle_buffer->processor_buffer_map.size());
    assert(maxCount == mc.particle_buffer->buffer_size/BUFFER_FRACTION);
    for (const auto &pair : mc.particle_buffer->processor_buffer_map) {
      assert(pair.second < nMessages);
      assert(ranks[pair.second] == pair.first);
    }
    return;
  }

  nMessages = mc.particle_buffer->num_buffers;
  maxCount = mc.particle_buffer->buffer_size/BUFFER_FRACTION;
  ranks = new int[nMessages];
  hipCalloc(nMessages,recvCounts);
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
  hipCalloc(nMessages*maxCount,recvParts);
  hipCalloc(nMessages,sendCounts);
  hipCalloc(nMessages*maxCount,sendParts);
  CHECK(hipEventCreate(&event));
}

void Messages::addParticle(const MC_Base_Particle &part, const int buffer)
{
  assert(buffer < nMessages);
  const int i = sendCounts[buffer]++;
  assert(i < maxCount);
  sendParts[buffer*maxCount+i] = part;
}

__global__ __launch_bounds__(1)
static void recvInit(const int nReceived, int *__restrict__ const nProcessing, int *__restrict__ const nExtras)
{
  *nProcessing = *nExtras+nReceived;
  *nExtras = 0;
}

__global__ __launch_bounds__(64)
static void recvCopy(const int nMessages, const int maxCount, const int totalCount, const int *const recvCounts, const MessageParticle *const recvParts, const int *const nParts, DeviceParticle *__restrict__ const parts)
{
  const int j = blockDim.x*blockIdx.x+threadIdx.x;
  int sumCount = 0;
  for (int i = 0; i < nMessages; i++) {
    const int count = recvCounts[i];
    sumCount += count;
    if (j < sumCount) {
      parts[*nParts-totalCount+j] = recvParts[i*maxCount+j-sumCount+count];
      return;
    }
  }
}

void Messages::completeRecvs(Device &device)
{
  MPI_Waitall(nMessages,recvReqs,recvStats);
  int total = 0;
  for (int i = 0; i < nMessages; i++) {
    int count = 0;
    MPI_Get_count(recvStats+i,mpiParticle,&count);
    recvCounts[i] = count;
    total += count;
  }
  std::swap(device.processing,device.extras);
  recvInit<<<1,1>>>(total,device.particleSizes+Device::PROCESSING,device.particleSizes+Device::EXTRAS);
  static constexpr int nt = 64;
  const int nb = (total+nt-1)/nt;
  recvCopy<<<nb,nt>>>(nMessages,maxCount,total,recvCounts,recvParts,device.particleSizes+Device::PROCESSING,device.processing);
}

void Messages::completeSends()
{
  for (int i = 0; i < nMessages; i++) sendCounts[i] = 0;
  MPI_Waitall(nMessages,sendReqs,MPI_STATUSES_IGNORE);
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
    const int count = sendCounts[i];
    assert(count <= maxCount);
    MPI_Isend(sendParts+offset,count,mpiParticle,ranks[i],tag,MPI_COMM_WORLD,sendReqs+i);
  }
}

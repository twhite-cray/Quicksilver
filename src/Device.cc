#include "Device.hh"
#include "MC_Base_Particle.hh"
#include "MC_Particle_Buffer.hh"
#include "MonteCarlo.hh"
#include "ParticleVaultContainer.hh"

Device::Device(MonteCarlo &mc): _mc(mc) {}

void Device::cycleInit()
{
  ParticleVaultContainer &pvc = *(_mc._particleVaultContainer);
  const MC_Particle_Buffer &pb = *(_mc.particle_buffer);

  const int vaultSize = pvc.getVaultSize();
  const int intSize = (MC_Base_Particle::num_base_ints * vaultSize + 2) * sizeof(int);
  const int doubleSize = MC_Base_Particle::num_base_floats * vaultSize * sizeof(double);
  const int charSize = MC_Base_Particle::num_base_chars * vaultSize * sizeof(char);

  _recvCount = intSize + doubleSize + charSize;

  const int nRecvs = pb.num_buffers;
  CHECK(hipHostMalloc(&_recvBuf, nRecvs * _recvCount));

  _recvSrcs.resize(nRecvs);
  _recvReqs.resize(nRecvs);

  for (int i = 0; i < nRecvs; i++) _recvSrcs[i] = pb.task[0].recv_buffer[i].processor;
}

void Device::cycleFinalize()
{
  for (MPI_Request &req : _recvReqs) MPI_Cancel(&req);
  CHECK(hipHostFree(_recvBuf));
  _recvBuf = NULL;
  _recvCount = 0;
}

static constexpr int tag = 33;

void Device::startRecvs()
{
  char *buf = _recvBuf;
  for (unsigned i = 0; i < _recvSrcs.size(); i++) {
    MPI_Irecv(buf, _recvCount, MPI_BYTE, _recvSrcs[i], tag, MPI_COMM_WORLD, _recvReqs.data()+i);
    buf += _recvCount;
  }
}

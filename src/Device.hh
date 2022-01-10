#pragma once
#include <mpi.h>
#include <vector>

class MonteCarlo;

class Device {
  public:
    Device(MonteCarlo &);
    void cycleInit();
    void cycleFinalize();
    void startRecvs();
  protected:
    MonteCarlo &_mc;

    char *_recvBuf;
    int _recvCount;
    std::vector<int> _recvSrcs;
    std::vector<MPI_Request> _recvReqs;
};

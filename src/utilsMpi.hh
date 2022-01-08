#ifndef UTILS_MPI_HH
#define UTILS_MPI_HH

#if defined (GNU_PERMISSIVE)
#pragma GCC diagnostic ignored "-fpermissive"
#endif

#include <mpi.h>

#if defined (GNU_PERMISSIVE)
#pragma GCC diagnostic ignored "-pedantic"
#endif

#ifndef MPI_INT64_T
#define MPI_INT64_T  MPI_LONG_LONG
#endif

#ifndef MPI_UINT64_T
#define MPI_UINT64_T MPI_UNSIGNED_LONG_LONG
#endif

double mpiWtime        ( void );
void mpiTest_cancelled ( MPI_Status *status, int *flag );
void mpiInit           ( int * argc, char *** argv );
void mpiFinalize       ( void );
void mpiComm_rank      ( MPI_Comm comm, int *rank );
void mpiComm_size      ( MPI_Comm comm, int *size );
int  mpiComm_split     ( MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
void mpiBarrier        ( MPI_Comm comm );
void mpiGet_version    ( int *version, int *subversion );
void mpiReduce         ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm );
void mpiGather         ( void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
void mpiBcast          ( void* buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
void mpiCancel         ( MPI_Request *request );
void mpiWait           ( MPI_Request *request, MPI_Status *status );
void mpiWaitall        ( int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses );
void mpiTest           ( MPI_Request *, int *, MPI_Status * );
void mpiIrecv          ( void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
void mpiRecv           ( void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
void mpiIsend          ( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
void mpiSend           ( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
void mpiType_contiguous( int count, MPI_Datatype old_type, MPI_Datatype *newtype );
void mpiType_commit    ( MPI_Datatype *datatype ) ;
void mpiAllreduce      ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op operation, MPI_Comm comm );
void mpiIAllreduce     ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op operation, MPI_Comm comm, MPI_Request *request);
void mpiScan           ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op operation, MPI_Comm comm );
void mpiAbort          ( MPI_Comm comm, int errorcode );
void mpiRequestFree    ( MPI_Request *request );

#endif  // end #ifndef UTILS_MPI_HH

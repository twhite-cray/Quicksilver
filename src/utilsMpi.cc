#include "utilsMpi.hh"
#include <cstdio>
#include <string.h>     // needed for memcpy on some compilers
#include <time.h>       // needed for clock
#include "qs_assert.hh"
#include "macros.hh"
#include "MonteCarlo.hh"
#include "MC_Processor_Info.hh"
#include "Globals.hh"


void mpiInit( int *argc, char ***argv)
{
   { // limit scope
      int err = MPI_Init(argc, argv);
      qs_assert(err == MPI_SUCCESS);
   } //limit scope
}


double mpiWtime( void ) { return MPI_Wtime(); }

int  mpiComm_split ( MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
   qs_assert(MPI_Comm_split(comm, color, key, newcomm) == MPI_SUCCESS); 
   return MPI_SUCCESS;
}

void mpiComm_rank( MPI_Comm comm, int *rank ) {   qs_assert(MPI_Comm_rank(comm, rank) == MPI_SUCCESS); }
void mpiCancel( MPI_Request *request ) { qs_assert(MPI_Cancel(request) == MPI_SUCCESS); }
void mpiTest_cancelled( MPI_Status *status, int *flag ) { qs_assert(MPI_Test_cancelled(status, flag) == MPI_SUCCESS); }
void mpiTest( MPI_Request *request, int *flag, MPI_Status * status) { qs_assert(MPI_Test(request, flag, status) == MPI_SUCCESS); }
void mpiWait( MPI_Request *request, MPI_Status *status ) { qs_assert(MPI_Wait(request, status) == MPI_SUCCESS); }
void mpiComm_size( MPI_Comm comm, int *size ) { qs_assert(MPI_Comm_size(comm, size) == MPI_SUCCESS); }
void mpiBarrier( MPI_Comm comm) { qs_assert(MPI_Barrier(comm) == MPI_SUCCESS); }
void mpiGet_version( int *version, int *subversion ) { qs_assert(MPI_Get_version(version, subversion) == MPI_SUCCESS); }
void mpiFinalize( void ) { qs_assert(MPI_Finalize() == MPI_SUCCESS); }
void mpiAbort( MPI_Comm comm, int errorcode ) { qs_assert(MPI_Abort(comm, errorcode) == MPI_SUCCESS); }
void mpiRequestFree( MPI_Request *request ){qs_assert( MPI_Request_free( request ) == MPI_SUCCESS);}

void mpiScan( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op operation, MPI_Comm comm )
   { qs_assert(MPI_Scan(sendbuf, recvbuf, count, datatype, operation, comm) == MPI_SUCCESS); }
void mpiType_commit(MPI_Datatype *datatype )
   { qs_assert(MPI_Type_commit( datatype ) == MPI_SUCCESS); }
void mpiType_contiguous(int count, MPI_Datatype old_type, MPI_Datatype *newtype)
   { qs_assert(MPI_Type_contiguous(count, old_type, newtype) == MPI_SUCCESS); }
void mpiWaitall( int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses )
   { qs_assert(MPI_Waitall(count, array_of_requests, array_of_statuses) == MPI_SUCCESS); }
void mpiAllreduce ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op operation, MPI_Comm comm )
   { qs_assert(MPI_Allreduce(sendbuf, recvbuf, count, datatype, operation, comm) == MPI_SUCCESS); }
void mpiIAllreduce( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op operation, MPI_Comm comm, MPI_Request *request)
   { qs_assert(MPI_Iallreduce(sendbuf, recvbuf, count, datatype, operation, comm, request) == MPI_SUCCESS); }
void mpiReduce( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm )
   { qs_assert(MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm) == MPI_SUCCESS); }
void mpiGather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
   { qs_assert(MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm) == MPI_SUCCESS); }
void mpiBcast( void* buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
   { qs_assert(MPI_Bcast(buf, count, datatype, root, comm) == MPI_SUCCESS); }
void mpiIrecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
   { qs_assert(MPI_Irecv(buf, count, datatype, source, tag, comm, request) == MPI_SUCCESS); }
void mpiRecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
   { qs_assert(MPI_Recv(buf, count, datatype, source, tag, comm, status) == MPI_SUCCESS); }
void mpiIsend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
   { qs_assert(MPI_Isend(buf, count, datatype, dest, tag, comm, request) == MPI_SUCCESS); }
void mpiSend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
   { qs_assert(MPI_Send(buf, count, datatype, dest, tag, comm) == MPI_SUCCESS); }
    


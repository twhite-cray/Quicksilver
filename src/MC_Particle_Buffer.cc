#include "MC_Particle_Buffer.hh"
#include <time.h>
#include "utilsMpi.hh"
#include "ParticleVaultContainer.hh"
#include "SendQueue.hh"
#include "MCT.hh"
#include "MC_Processor_Info.hh"
#include "Tallies.hh"
#include "MonteCarlo.hh"
#include "Globals.hh"
#include "MC_Fast_Timer.hh"
#include "macros.hh"
#include "NVTX_Range.hh"

//
//  mcp_test_done_class
//

//----------------------------------------------------------------------------------------------------------------------
//  Reset to 0 and clear all vectors.
//----------------------------------------------------------------------------------------------------------------------
void mcp_test_done_class::Zero_Out()
{
    this->local_sent = 0;
    this->local_recv = 0;
    this->BlockingSum = 0;
}

//----------------------------------------------------------------------------------------------------------------------
//  Free the memory and cleanup communication.
//----------------------------------------------------------------------------------------------------------------------
void mcp_test_done_class::Free_Memory()
{
    this->Zero_Out();
}

//----------------------------------------------------------------------------------------------------------------------
//  Get the number of particles created and completed.
//----------------------------------------------------------------------------------------------------------------------
void mcp_test_done_class::Get_Local_Gains_And_Losses(MonteCarlo *monteCarlo, int64_t sent_recv[2])
{
    uint64_t gains = 0, losses = 0;

    Balance &bal = monteCarlo->_tallies->_balanceTask[0]; // SumTasks has been called, so just use index 0

    gains   = bal._start  + bal._source + bal._produce + bal._split;
    losses  = bal._absorb + bal._census + bal._escape  + bal._rr;
    losses += bal._fission; 

    sent_recv[0] = gains;
    sent_recv[1] = losses;
}

//
// MC_Particle_Buffer :: PRIVATE Functions
//

//----------------------------------------------------------------------------------------------------------------------
//  Initializes the particle buffers, mallocing them and assigning processors to buffers.
//
//----------------------------------------------------------------------------------------------------------------------
void MC_Particle_Buffer::Instantiate()
{
    this->test_done.Zero_Out();
    this->Initialize_Map();
}

//----------------------------------------------------------------------------------------------------------------------
//  Define this->processor_buffer_map[...].
//
//----------------------------------------------------------------------------------------------------------------------
void MC_Particle_Buffer::Initialize_Map()
{
    this->num_buffers = 0;

    // Determine number of buffers needed and assign processors to buffers
    for ( int domain_index = 0; domain_index < mcco->domain.size(); domain_index++ )
    {
        MC_Domain &domain = mcco->domain[domain_index];
        for ( int neighbor_index = 0; neighbor_index < domain.mesh._nbrRank.size(); neighbor_index++ )
        {
            int neighbor_rank = domain.mesh._nbrRank[neighbor_index];

            // If neighbor is not on same processor
            if ( neighbor_rank != mcco->processor_info->rank )
            {
                if ( this->Choose_Buffer(neighbor_rank) == -1 )
                {
                    this->processor_buffer_map[neighbor_rank] = this->num_buffers++;
                }
            }
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------
//  Trivially_Done
//      Do we have more than 1 processor running
//      Do we have any particles to process locally 
//----------------------------------------------------------------------------------------------------------------------
bool MC_Particle_Buffer::Trivially_Done()
{
    if (mcco->processor_info->num_processors > 1) 
    {
        return false;
    }

    uint64_t processingSize = mcco->_particleVaultContainer->sizeProcessing();
    if( processingSize == 0 )
    {
        return true;
    }

    return false;
}

//
// MC_Particle_Buffer :: PUBLIC Functions
// 

//----------------------------------------------------------------------------------------------------------------------
//  Constructor.
//----------------------------------------------------------------------------------------------------------------------
MC_Particle_Buffer::MC_Particle_Buffer(MonteCarlo *mcco_, size_t bufferSize_)
{
    this->mcco  = mcco_;
    this->test_done.Zero_Out();
    this->num_buffers = 0;
    this->buffer_size = bufferSize_;
    this->processor_buffer_map.clear();
}

//----------------------------------------------------------------------------------------------------------------------
//  Initializes the particle buffers, mallocing them and assigning processors to buffers.
//
//----------------------------------------------------------------------------------------------------------------------
void MC_Particle_Buffer::Initialize()
{
    NVTX_Range range("MC_Particle_Buffer::Initialize");

    if (mcco->processor_info->num_processors > 1) 
    {
        this->Instantiate();

        mpiBarrier(mcco->processor_info->comm_mc_world);
    }
}

//----------------------------------------------------------------------------------------------------------------------
//  Selects a  particle buffer.
//----------------------------------------------------------------------------------------------------------------------
int MC_Particle_Buffer::Choose_Buffer(int processor)
{
    std::map<int,int>::iterator it = this->processor_buffer_map.find(processor);

    if ( it == this->processor_buffer_map.end() )
    {
        // return -1 if the input processor does not communicate with this processor. i.e. not found in map.
        return -1;
    }
    else
    {
        return (*it).second;
    }
}

//----------------------------------------------------------------------------------------------------------------------
//  Test to see if we are done with streaming communication.
//----------------------------------------------------------------------------------------------------------------------
bool MC_Particle_Buffer::Test_Done_New() 
{
    if ( !(mcco->processor_info->num_processors > 1 ))
    {
        return this->Trivially_Done();
    }

    MC_VERIFY_THREAD_ZERO

    MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_Test_Done);

    mcco->_tallies->SumTasks();

    bool answer = this->Allreduce_ParticleCounts();

    MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_Test_Done);
    return answer;
}

//----------------------------------------------------------------------------------------------------------------------
//  Perform blocking allreduce on particle counts, return true if counts are equal.
//----------------------------------------------------------------------------------------------------------------------
bool MC_Particle_Buffer::Allreduce_ParticleCounts()
{
    int64_t buf[2];
    int64_t hard_blocking_sum[2] = {0, 0};
    this->test_done.Get_Local_Gains_And_Losses(mcco, buf);
    mpiAllreduce(buf, hard_blocking_sum, 2, MPI_INT64_T, MPI_SUM, mcco->processor_info->comm_mc_world);
    return ( hard_blocking_sum[0] == hard_blocking_sum[1] );
}


//----------------------------------------------------------------------------------------------------------------------
//  Free the particle buffer memory.
//----------------------------------------------------------------------------------------------------------------------
void MC_Particle_Buffer::Free_Memory()
{
    this->num_buffers = 0;      // buffers are now freed
    this->processor_buffer_map.clear();
    this->test_done.Free_Memory();
}


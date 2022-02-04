#include <iostream>
#include "utils.hh"
#include "Parameters.hh"
#include "utilsMpi.hh"
#include "MonteCarlo.hh"
#include "initMC.hh"
#include "Tallies.hh"
#include "PopulationControl.hh"
#include "ParticleVaultContainer.hh"
#include "ParticleVault.hh"
#include "MC_Particle_Buffer.hh"
#include "MC_Processor_Info.hh"
#include "MC_Time_Info.hh"
#include "macros.hh"
#include "MC_Fast_Timer.hh"
#include "MC_SourceNow.hh"
#include "SendQueue.hh"
#include "NVTX_Range.hh"
#include "cudaUtils.hh"
#include "cudaFunctions.hh"
#include "qs_assert.hh"
#include "CycleTracking.hh"
#include "CoralBenchmark.hh"
#include "EnergySpectrum.hh"

#include "git_hash.hh"
#include "git_vers.hh"

void gameOver();
void cycleInit( bool loadBalance );
void cycleTracking(MonteCarlo* monteCarlo);
void cycleFinalize();

using namespace std;

MonteCarlo *mcco  = NULL;

int main(int argc, char** argv)
{
   mpiInit(&argc, &argv);
   printBanner(GIT_VERS, GIT_HASH);

   Parameters params = getParameters(argc, argv);
   printParameters(params, cout);

   // mcco stores just about everything. 
   mcco = initMC(params); 

   int loadBalance = params.simulationParams.loadBalance;

   MC_FASTTIMER_START(MC_Fast_Timer::main);     // this can be done once mcco exist.

   const int nSteps = params.simulationParams.nSteps;

   for (int ii=0; ii<nSteps; ++ii)
   {
      cycleInit( bool(loadBalance) );
      cycleTracking(mcco);
      cycleFinalize();

      mcco->fast_timer->Last_Cycle_Report(
            params.simulationParams.cycleTimers,
            mcco->processor_info->rank,
            mcco->processor_info->num_processors,
            mcco->processor_info->comm_mc_world );
   }


   MC_FASTTIMER_STOP(MC_Fast_Timer::main);

   gameOver();

   coralBenchmarkCorrectness(mcco, params);

   delete mcco;

   mpiFinalize();
   
   return 0;
}

void gameOver()
{
    mcco->fast_timer->Cumulative_Report(mcco->processor_info->rank,
                                        mcco->processor_info-> num_processors,
                                        mcco->processor_info->comm_mc_world,
                                        mcco->_tallies->_balanceCumulative._numSegments);
    mcco->_tallies->_spectrum.PrintSpectrum(mcco);
}

void cycleInit( bool loadBalance )
{

    MC_FASTTIMER_START(MC_Fast_Timer::cycleInit);

    mcco->clearCrossSectionCache();

    mcco->_tallies->CycleInitialize(mcco);

    mcco->_particleVaultContainer->swapProcessingProcessedVaults();

    mcco->_tallies->_balanceTask[0]._start = mcco->_particleVaultContainer->sizeProcessing();

    mcco->particle_buffer->Initialize();
    mcco->_messages.init();

    MC_SourceNow(mcco);
   
    PopulationControl(mcco, loadBalance); // controls particle population

    RouletteLowWeightParticles(mcco); // Delete particles with low statistical weight

    mcco->_device.cycleInit(*mcco);

    MC_FASTTIMER_STOP(MC_Fast_Timer::cycleInit);
}


void cycleTracking(MonteCarlo *monteCarlo)
{
    MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking);

    bool done = false;
    Device &device = monteCarlo->_device;
    Messages &messages = monteCarlo->_messages;

    ParticleVaultContainer &my_particle_vault = *(monteCarlo->_particleVaultContainer);

    do
    {
        int particle_count = 0; // Initialize count of num_particles processed

        while ( !done )
        {
            messages.startRecvs();
            MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_Kernel);

            ParticleVault *processingVault = my_particle_vault.getTaskProcessingVault();
            ParticleVault *processedVault =  my_particle_vault.getTaskProcessedVault();

            int numParticles = processingVault->size();
            *device.processingSize = numParticles;
            for (int i = 0; i < numParticles; i++) device.processing[i] = (*processingVault)[i];

            if ( numParticles != 0 )
            {
              NVTX_Range trackingKernel("cycleTracking_TrackingKernel"); 
              CycleTrackingGuts( monteCarlo, numParticles, processingVault, processedVault, device );
            }

            for (int i = 0; i < numParticles; i++) assert(device.processing[i] == (*processingVault)[i]);

            particle_count += numParticles;

            MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_Kernel);

            MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_MPI);

            NVTX_Range cleanAndComm("cycleTracking_clean_and_comm");

            SendQueue &sendQueue = *(my_particle_vault.getSendQueue());
            for ( int index = 0; index < sendQueue.size(); index++ )
            {
              sendQueueTuple& sendQueueT = sendQueue.getTuple( index );
              MC_Base_Particle mcb_particle;

              processingVault->getBaseParticleComm( mcb_particle, sendQueueT._particleIndex );

              int buffer = monteCarlo->particle_buffer->Choose_Buffer(sendQueueT._neighbor );
              messages.addParticle(mcb_particle,buffer);
            }

            messages.startSends();
            processingVault->clear(); //remove the invalid particles
            sendQueue.clear();

            // Move particles in "extra" vault into the regular vaults.
            my_particle_vault.cleanExtraVault();
            MPI_Barrier(MPI_COMM_WORLD);

            // receive any particles that have arrived from other ranks
            messages.completeRecvs();
            messages.completeSends();

            MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_MPI);

            MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_MPI);

            //Test for done - blocking on all MPI ranks
            NVTX_Range doneRange("cycleTracking_Test_Done_New");
            done = monteCarlo->particle_buffer->Test_Done_New();
            doneRange.endRange();

            MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_MPI);

        } // while not done: Test_Done_New()

        // Everything should be done normally.
        done = monteCarlo->particle_buffer->Test_Done_New();

    } while ( !done );

   MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking);
}


void cycleFinalize()
{
    MC_FASTTIMER_START(MC_Fast_Timer::cycleFinalize);

    mcco->_device.cycleFinalize(*mcco);

    mcco->_tallies->_balanceTask[0]._end = mcco->_particleVaultContainer->sizeProcessed();

    // Update the cumulative tally data.
    mcco->_tallies->CycleFinalize(mcco); 

    mcco->time_info->cycle++;

    mcco->particle_buffer->Free_Memory();

    MC_FASTTIMER_STOP(MC_Fast_Timer::cycleFinalize);
}


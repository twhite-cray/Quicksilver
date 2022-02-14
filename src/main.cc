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
    mcco->_messagesA.init(*mcco);
    mcco->_messagesB.init(*mcco);

    MC_SourceNow(mcco);
   
    PopulationControl(mcco, loadBalance); // controls particle population

    RouletteLowWeightParticles(mcco); // Delete particles with low statistical weight

    mcco->_device.cycleInit(*mcco);

    MC_FASTTIMER_STOP(MC_Fast_Timer::cycleInit);
}


void cycleTracking(MonteCarlo *monteCarlo)
{
    MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking);

    Device &device = monteCarlo->_device;
    Messages &ma = monteCarlo->_messagesA;
    Messages &mb = monteCarlo->_messagesB;

    const int nMid = device.particleSizes[Device::PROCESSING]/2;
    ma.startRecvs();
    CycleTrackingGuts(0,nMid,device,ma.maxCount,ma.sendCounts,ma.sendParts);
    mb.startRecvs();
    CycleTrackingGuts(nMid,-1,device,mb.maxCount,mb.sendCounts,mb.sendParts);
    bool doA = true;
    bool doB = true;
    while (doA || doB) {
      if (doA) {
        ma.startSends();
        ma.completeRecvs(device);
        ma.completeSends();
        doA = !monteCarlo->particle_buffer->Test_Done_New();
      }
      if (doA) {
        ma.startRecvs();
        CycleTrackingGuts(0,-1,device,ma.maxCount,ma.sendCounts,ma.sendParts);
      }
      if (doB) {
        mb.startSends();
        mb.completeRecvs(device);
        mb.completeSends();
        doB = !monteCarlo->particle_buffer->Test_Done_New();
      }
      if (doB) {
        mb.startRecvs();
        CycleTrackingGuts(0,-1,device,mb.maxCount,mb.sendCounts,mb.sendParts);
      }
    }

    MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking);
}


void cycleFinalize()
{
    MC_FASTTIMER_START(MC_Fast_Timer::cycleFinalize);

    mcco->_device.cycleFinalize(*mcco);

    assert(mcco->_particleVaultContainer->sizeProcessed() == mcco->_device.particleSizes[Device::PROCESSED]);
    mcco->_tallies->_balanceTask[0]._end = mcco->_particleVaultContainer->sizeProcessed();

    // Update the cumulative tally data.
    mcco->_tallies->CycleFinalize(mcco); 

    mcco->time_info->cycle++;

    mcco->particle_buffer->Free_Memory();

    MC_FASTTIMER_STOP(MC_Fast_Timer::cycleFinalize);
}


#include <iostream>
#include "utils.hh"
#include "Parameters.hh"
#include "utilsMpi.hh"
#include "CollisionEvent.hh"
#include "MonteCarlo.hh"
#include "initMC.hh"
#include "Tallies.hh"
#include "PopulationControl.hh"
#include "ParticleVaultContainer.hh"
#include "ParticleVault.hh"
#include "MC_Facet_Crossing_Event.hh"
#include "MC_Particle_Buffer.hh"
#include "MC_Processor_Info.hh"
#include "MC_Segment_Outcome.hh"
#include "MC_Time_Info.hh"
#include "macros.hh"
#include "MC_Fast_Timer.hh"
#include "MC_SourceNow.hh"
#include "MCT.hh"
#include "NVTX_Range.hh"
#include "cudaUtils.hh"
#include "cudaFunctions.hh"
#include "qs_assert.hh"
#include "CoralBenchmark.hh"
#include "EnergySpectrum.hh"

#include "git_hash.hh"
#include "git_vers.hh"

static constexpr int NT = 64;

__global__ __launch_bounds__(NT) static void CycleTrackingGuts( const int ipLo, int ipHi, Device device, const int maxCount, int *__restrict__ const sendCounts, MessageParticle *__restrict__ const sendParts)
{
    const int di = gridDim.x*blockDim.x;
    ipHi = (ipHi < 0) ? device.particleSizes[Device::PROCESSING] : ipHi;
    MC_Particle mc_particle;
    int ipOld = ipLo-1;
    int ip = ipLo+blockIdx.x*blockDim.x+threadIdx.x;
    while (ip < ipHi) {

        if (ipOld != ip) {
          ipOld = ip;
          device.processing[ip].set(mc_particle);
          if (mc_particle.time_to_census <= 0) mc_particle.time_to_census += device.timeStep;
          mc_particle.energy_group = device.getEnergyGroup(mc_particle.kinetic_energy);
        }

        // Determine the outcome of a particle at the end of this segment such as:
        //
        //   (0) Undergo a collision within the current cell,
        //   (1) Cross a facet of the current cell,
        //   (2) Reach the end of the time step and enter census,
        //
        const MC_Segment_Outcome_type::Enum segment_outcome = MC_Segment_Outcome(device, mc_particle);

        atomicFetchAdd(device.tallies+Device::SEGMENTS,1UL);

        mc_particle.num_segments += 1.;  /* Track the number of segments this particle has
                                            undergone this cycle on all processes. */
        switch (segment_outcome) {
        case MC_Segment_Outcome_type::Collision:
            {
              if (CollisionEvent(device, mc_particle) != MC_Collision_Event_Return::Continue_Tracking) ip += di;
            }
            break;
        case MC_Segment_Outcome_type::Facet_Crossing:
            {
                // The particle has reached a cell facet.
                const MC_Tally_Event::Enum facet_crossing_type = MC_Facet_Crossing_Event(mc_particle, device, ip, maxCount, sendCounts, sendParts);
                if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Transit_Exit)
                {}
                else if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Escape)
                {
                    atomicFetchAdd(device.tallies+Device::ESCAPE,1UL);
                    mc_particle.last_event = MC_Tally_Event::Facet_Crossing_Escape;
                    ip += di;
                }
                else if (facet_crossing_type == MC_Tally_Event::Facet_Crossing_Reflection)
                {
                    MCT_Reflect_Particle(device, mc_particle);
                }
                else
                {
                    // Enters an adjacent cell in an off-processor domain.
                    ip += di;
                }
            }
            break;
    
        case MC_Segment_Outcome_type::Census:
            {
                const int iProcessed = atomicFetchAdd(device.particleSizes+Device::PROCESSED,1);
                device.processed[iProcessed] = mc_particle;
                atomicFetchAdd(device.tallies+Device::CENSUS,1UL);
                ip += di;
                break;
            }
            
        default:
           abort();
           break;  // should this be an error
        }
    }
}
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

    //const int nb = (monteCarlo->processor_info->thread_target+NT-1)/NT;

    Device &device = monteCarlo->_device;
    Messages &ma = monteCarlo->_messagesA;
    Messages &mb = monteCarlo->_messagesB;

    const int nb = (monteCarlo->processor_info->thread_target+NT-1)/NT;

    const int nMid = device.particleSizes[Device::PROCESSING]/2;
    ma.startRecvs();
    CycleTrackingGuts<<<nb,NT>>>(0,nMid,device,ma.maxCount,ma.sendCounts,ma.sendParts);
    CHECK(hipDeviceSynchronize());
    mb.startRecvs();
    CycleTrackingGuts<<<nb,NT>>>(nMid,-1,device,mb.maxCount,mb.sendCounts,mb.sendParts);
    CHECK(hipDeviceSynchronize());
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
        CycleTrackingGuts<<<nb,NT>>>(0,-1,device,ma.maxCount,ma.sendCounts,ma.sendParts);
        CHECK(hipDeviceSynchronize());
      }
      if (doB) {
        mb.startSends();
        mb.completeRecvs(device);
        mb.completeSends();
        doB = !monteCarlo->particle_buffer->Test_Done_New();
      }
      if (doB) {
        mb.startRecvs();
        CycleTrackingGuts<<<nb,NT>>>(0,-1,device,mb.maxCount,mb.sendCounts,mb.sendParts);
        CHECK(hipDeviceSynchronize());
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


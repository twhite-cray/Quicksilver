#include "MonteCarlo.hh"
#include "NuclearData.hh"
#include "MaterialDatabase.hh"
#include "ParticleVaultContainer.hh"
#include "MC_RNG_State.hh"
#include "Tallies.hh"
#include "MC_Processor_Info.hh"
#include "MC_Time_Info.hh"
#include "MC_Particle_Buffer.hh"
#include "MC_Fast_Timer.hh"
#include <cmath>

#include "macros.hh" // current location of openMP wrappers.
#include "cudaUtils.hh"

using std::ceil;

//----------------------------------------------------------------------------------------------------------------------
// Construct a MonteCarlo object.
//----------------------------------------------------------------------------------------------------------------------
MonteCarlo::MonteCarlo(const Parameters& params)
: _params(params),
  _nuclearData(NULL)
{
   _nuclearData            = 0;
   _materialDatabase       = 0;

   _tallies                = new Tallies( params.simulationParams.balanceTallyReplications, 
       params.simulationParams.fluxTallyReplications,
       params.simulationParams.cellTallyReplications,
       params.simulationParams.energySpectrum,
       params.simulationParams.nGroups);
   processor_info          = new MC_Processor_Info();
   time_info               = new MC_Time_Info();
   fast_timer              = new MC_Fast_Timer_Container();

   source_particle_weight = 0.0;

    size_t num_processors = processor_info->num_processors;
    size_t num_particles  = params.simulationParams.nParticles;

    size_t num_particles_on_process = num_particles / num_processors;

    if( num_particles_on_process <= 0 )
    {
        MC_Fatal_Jump( "Not enough particles for each process ( Ranks: %d Num Particles: %d ) \n", num_processors, num_particles ); 
        num_particles_on_process = 1;
    }

    size_t vector_size = 0;

    for (auto matIter  = params.materialParams.begin(); 
              matIter != params.materialParams.end(); 
              matIter++)
    {
        const MaterialParameters& mp = matIter->second;
        double nuBar = params.crossSectionParams.at(mp.fissionCrossSection).nuBar;
        size_t nb = ceil( nuBar );
        size_t test_size = nb*( num_particles_on_process );

        if ( test_size > vector_size )
            vector_size = test_size;
    }
    if ( vector_size == 0 )
        vector_size = 2*num_particles_on_process;

    int num_extra_vaults = ( vector_size / num_particles_on_process ) + 1;
    //Previous definition was not enough extra space for some reason? need to determine why still

    particle_buffer         = new MC_Particle_Buffer(this, num_particles_on_process);
    _particleVaultContainer = new ParticleVaultContainer(num_particles_on_process, 1, num_extra_vaults);
}

//----------------------------------------------------------------------------------------------------------------------
// Destruct a MonteCarlo object.
//----------------------------------------------------------------------------------------------------------------------
MonteCarlo::~MonteCarlo()
{
  delete _nuclearData;
  delete _particleVaultContainer;
  delete _materialDatabase;
  delete _tallies;
  delete processor_info;
  delete time_info;
  delete fast_timer;
  delete particle_buffer;
}

void MonteCarlo::clearCrossSectionCache()
{
   int numEnergyGroups = _nuclearData->_numEnergyGroups;
   for (unsigned ii=0; ii<domain.size(); ++ii)
      domain[ii].clearCrossSectionCache(numEnergyGroups);
}


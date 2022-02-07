#ifndef PARTICLEVAULTCONTAINER_HH
#define PARTICLEVAULTCONTAINER_HH

#include "DeclareMacro.hh"

#include "portability.hh"
#include "ParticleVault.hh"
#include "SendQueue.hh"

//---------------------------------------------------------------
// ParticleVaultContainer is a container of ParticleVaults. 
// These Vaults are broken down into user defined chunks that can 
// be used to overlap asynchronous MPI with the tracking kernel.
//
// Facilities for storing Processing, Processed, and Extra vaults 
// are controled by the ParticleVaultContainer. As well as the 
// sendQueue, which lists the particles that must be send to 
// another process via MPI
//--------------------------------------------------------------

class MC_Base_Particle;
class MC_Particle;

typedef unsigned long long int uint64_cu;

class ParticleVaultContainer
{
  public:
    
    //Constructor
    ParticleVaultContainer( uint64_t vault_size );

    //Destructor
    ~ParticleVaultContainer();

    //Basic Getters
    uint64_t getVaultSize(){      return _vaultSize; }

    ParticleVault *getTaskProcessingVault() { return &_processingVault; }
    ParticleVault *getTaskProcessedVault() { return &_processedVault; }

    //Returns a pointer to the Send Queue
    SendQueue* getSendQueue();

    //Counts Particles in all vaults
    uint64_t sizeProcessing() const { return _processingVault.size(); }
    uint64_t sizeProcessed() const { return _processedVault.size(); }

    //Swaps the particles in Processed for the empty vaults in 
    //Processing
    void swapProcessingProcessedVaults();

    //Adds a particle to the processing particle vault
    void addProcessingParticle( MC_Base_Particle &particle );
    //Adds a particle to the extra particle vault
    void addExtraParticle( MC_Particle &particle );
 
    //Pushes particles from Extra Vault onto the Processing 
    //Vault list
    void cleanExtraVault();

  private:
    
    //The Size of the ParticleVaults (fixed at runtime for 
    //each run)
    uint64_t _vaultSize;

    //The send queue - stores particle index and neighbor index 
    //for any particles that hit (TRANSIT_OFF_PROCESSOR) 
    SendQueue _sendQueue;

  public:

    ParticleVault _processingVault;
    ParticleVault _processedVault;
    ParticleVault _extraVault;
};

#endif

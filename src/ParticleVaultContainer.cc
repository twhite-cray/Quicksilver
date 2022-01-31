#include "ParticleVaultContainer.hh"
#include "ParticleVault.hh"
#include "SendQueue.hh"
#include "MemoryControl.hh"
#include "qs_assert.hh"

//--------------------------------------------------------------
//------------ParticleVaultContainer Constructor----------------
//Sets up the fixed sized data and pre-allocates the minimum 
//needed for processing and processed vaults
//--------------------------------------------------------------

ParticleVaultContainer::
ParticleVaultContainer( uint64_t vault_size )
: _vaultSize      ( vault_size       )
{
    _processedVault.reserve( vault_size );
    _processingVault.reserve( vault_size );
    _extraVault.reserve( vault_size );

    _sendQueue = MemoryControl::allocate<SendQueue>(1 ,VAR_MEM);
    _sendQueue->reserve( vault_size );
}

//--------------------------------------------------------------
//------------ParticleVaultContainer Destructor-----------------
//Deletes memory allocaetd using the Memory Control class
//--------------------------------------------------------------

ParticleVaultContainer::
~ParticleVaultContainer()
{
    MemoryControl::deallocate( _sendQueue, 1, VAR_MEM );
}

//--------------------------------------------------------------
//------------getSendQueue--------------------------------------
//Returns a pointer to the Send Queue
//--------------------------------------------------------------
HOST_DEVICE
SendQueue* ParticleVaultContainer::
getSendQueue()
{
    return this->_sendQueue;
}
HOST_DEVICE_END

//--------------------------------------------------------------
//------------swapProcessingProcessedVaults---------------------
//Swaps the vaults from Processed that have particles in them
//with empty vaults from processing to prepair for the next
//cycle
//
//ASSUMPTIONS:: 
//  2) _processingVault is always empty of particles when this is
//      called
//--------------------------------------------------------------

void ParticleVaultContainer::
swapProcessingProcessedVaults()
{
    if (this->_processedVault.size() > 0) {
      std::swap(this->_processingVault, this->_processedVault);
    }
}

//--------------------------------------------------------------
//------------addProcessingParticle-----------------------------
//Adds a particle to the processing particle vault
//--------------------------------------------------------------

void ParticleVaultContainer::
addProcessingParticle( MC_Base_Particle &particle )
{
    _processingVault.pushBaseParticle(particle);
}

//--------------------------------------------------------------
//------------addExtraParticle----------------------------------
//adds a particle to the extra particle vaults (used in kernel)
//--------------------------------------------------------------
HOST_DEVICE
void ParticleVaultContainer::
addExtraParticle( MC_Particle &particle)
{
    _extraVault.pushParticle( particle );
}
HOST_DEVICE_END

//--------------------------------------------------------------
//------------cleanExtraVault----------------------------------
//Moves the particles from the _extraVault into the 
//_processedVault
//--------------------------------------------------------------

void ParticleVaultContainer::
cleanExtraVault()
{
  uint64_t size_extra = this->_extraVault.size();
  if( size_extra > 0 )
  {
    const uint64_t fill_size = this->_vaultSize - this->_processingVault.size();
    assert(size_extra < fill_size);
    this->_processingVault.collapse( fill_size, &(this->_extraVault) );
  }
}


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
ParticleVaultContainer( uint64_t vault_size, 
                        uint64_t num_vaults )
: _vaultSize      ( vault_size       ),
  _processedVault ( num_vaults       )
{
    //Allocate and reserve space for particles for each vault
    for( uint64_t vault = 0; vault < num_vaults; vault++ )
    {
        //Allocate Processed Vault
        _processedVault[vault]  = 
            MemoryControl::allocate<ParticleVault>(1 ,VAR_MEM);
        _processedVault[vault]->reserve( vault_size );
    }

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
    for( int64_t jj = _processedVault.size()-1; jj >= 0; jj-- )
    {
        MemoryControl::deallocate(_processedVault[jj], 1, VAR_MEM);
    }
    MemoryControl::deallocate( _sendQueue, 1, VAR_MEM );
}

//--------------------------------------------------------------
//------------getTaskProcessedVault-----------------------------
//Returns a pointer to the Particle Vault in the processed list
//at the index provided
//--------------------------------------------------------------

ParticleVault* ParticleVaultContainer::
getTaskProcessedVault(uint64_t vaultIndex)
{
//   qs_assert(vaultIndex >= 0);
//   qs_assert(vaultIndex < _processedVault.size());
   return _processedVault[vaultIndex];
}

//--------------------------------------------------------------
//------------getFirstEmptyProcessedVault-----------------------
//Returns a pointer to the first empty Particle Vault in the 
//processed list
//--------------------------------------------------------------

uint64_t ParticleVaultContainer::
getFirstEmptyProcessedVault()
{
    return 0;
    uint64_t index = 0;

    while( _processedVault[index]->size() != 0 )
    {
        index++;
        if( index == _processedVault.size() )
        {
            ParticleVault* vault = MemoryControl::allocate<ParticleVault>(1,VAR_MEM);
            vault->reserve( _vaultSize );
            this->_processedVault.push_back(vault);
            printf("TREY getFirstEmptyProcessedVault _processedVault %lu %lu\n",this->_processedVault.size(),index);
        }
    }
    return index;
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
//------------sizeProcessed-------------------------------------
//returns the total number of particles in the processed vault
//--------------------------------------------------------------

uint64_t ParticleVaultContainer::
sizeProcessed()
{
    uint64_t sum_size = 0;
    for( uint64_t vault = 0; vault < _processedVault.size(); vault++ )
    {
        sum_size += _processedVault[vault]->size();
    }
    return sum_size;
}

//--------------------------------------------------------------
//------------collapseProcessed---------------------------------
//Collapses the particles in the processed vault down to the
//first vaults needed to hold that number of particles
//--------------------------------------------------------------
    
void ParticleVaultContainer::
collapseProcessed()
{
    uint64_t num_vaults = this->_processedVault.size();

    uint64_t fill_vault_index = 0;
    uint64_t from_vault_index = num_vaults-1;

    while( fill_vault_index < from_vault_index )
    {
        if( _processedVault[fill_vault_index]->size() == this->_vaultSize )
        {
            fill_vault_index++;
        }
        else
        {
            if( this->_processedVault[from_vault_index]->size() == 0 )
            {
                from_vault_index--;
            }
            else
            {
                uint64_t fill_size = this->_vaultSize - this->_processedVault[fill_vault_index]->size();

                this->_processedVault[fill_vault_index]->collapse( fill_size, this->_processedVault[from_vault_index] );
            }
        }
    }
}

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
    assert(this->_processedVault.size() == 1);
    if (this->_processedVault[0]->size() > 0) {
      std::swap(this->_processingVault, *(this->_processedVault[0]));
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


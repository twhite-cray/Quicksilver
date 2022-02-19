#include "Device.hh"

#include "cudaUtils.hh"
#include "MC_Base_Particle.hh"
#include "MC_Particle_Buffer.hh"
#include "MC_Time_Info.hh"
#include "MonteCarlo.hh"
#include "NuclearData.hh"
#include "ParticleVaultContainer.hh"
#include "Tallies.hh"

void Device::init(MonteCarlo &mc)
{
  assert(domains == nullptr);
  const int domainSize = mc.domain.size();
  hipCalloc(domainSize,domains);

  const int groupSize = mc._nuclearData->_numEnergyGroups;

  int cellSizeSum = 0;
  int nodeSizeSum = 0;
  for (int i = 0; i < domainSize; i++) {
    cellSizeSum += mc.domain[i].cell_state.size();
    nodeSizeSum += mc.domain[i].mesh._node.size();
  }
  
  DeviceCell *cells = nullptr;
  hipCalloc(cellSizeSum,cells);
  double *totals = nullptr;
  hipCalloc(cellSizeSum*groupSize,totals);
  double *groupTallies = nullptr;
  hipCalloc(cellSizeSum*groupSize,groupTallies);
  double3 *nodes = nullptr;
  hipCalloc(nodeSizeSum,nodes);
  const int neighborSize = mc.domain[0].mesh._nbrRank.size();
  int *neighbors = nullptr;
  hipCalloc(neighborSize*domainSize,neighbors);
  for (int i = 0; i < domainSize; i++) {
    domains[i].cells = cells;
    const int cellSize = mc.domain[i].cell_state.size();
    for (int j = 0; j < cellSize; j++) {
      cells[j] = mc.domain[i].cell_state[j];
      cells[j].totals = totals;
      totals += groupSize;
      cells[j].groupTallies = groupTallies;
      groupTallies += groupSize;
      assert(DeviceCell::numFacets == mc.domain[i].mesh._cellConnectivity[j].num_facets);
      for (int k = 0; k < DeviceCell::numFacets; k++) {
        const MC_General_Plane &plane = mc.domain[i].mesh._cellGeometry[j]._facet[k];
        DeviceFacet &facet = cells[j].facets[k];
        facet.plane = double4{plane.A,plane.B,plane.C,plane.D};
        const int *const p = mc.domain[i].mesh._cellConnectivity[j]._facet[k].point;
        facet.point = int3{p[0],p[1],p[2]};
        const Subfacet_Adjacency &adjacency = mc.domain[i].mesh._cellConnectivity[j]._facet[k].subfacet;
        facet.adjacentCell = adjacency.adjacent.cell;
        facet.adjacentDomain = adjacency.adjacent.domain;
        facet.adjacentFacet = adjacency.adjacent.facet;
        facet.currentDomain = adjacency.current.domain;
        facet.neighbor = adjacency.neighbor_index;
        facet.event = adjacency.event;
      }
      for (int k = 0; k < DeviceCell::numQuadPoints; k++) {
        cells[j].quadPoints[k] = mc.domain[i].mesh._cellConnectivity[j]._point[k];
      }
    }
    cells += cellSize;
    domains[i].nodes = nodes;
    const int nodeSize = mc.domain[i].mesh._node.size();
    for (int j = 0; j < nodeSize; j++) {
      const MC_Vector &node = mc.domain[i].mesh._node[j];
      nodes[j] = double3{node.x,node.y,node.z};
    }
    nodes += nodeSize;
    domains[i].neighbors = neighbors;
    neighbors += neighborSize;
    assert(mc.domain[i].mesh._nbrRank.size() == neighborSize);
  }

  assert(mats == nullptr);
  const int matSize = mc._materialDatabase->_mat.size();
  hipCalloc(matSize,mats);
  
  int isoSizeSum = 0;
  for (int i = 0; i < matSize; i++) isoSizeSum += mc._materialDatabase->_mat[i]._iso.size();

  DeviceIsotope *isos = nullptr;
  hipCalloc(isoSizeSum,isos);
  for (int i = 0; i < matSize; i++) {
    mats[i].isos = isos;
    const int isoSize = mc._materialDatabase->_mat[i]._iso.size();
    mats[i].isoSize = isoSize;
    for (int j = 0; j < isoSize; j++) isos[j] = mc._materialDatabase->_mat[i]._iso[j];
    isos += isoSize;
    mats[i].mass = mc._materialDatabase->_mat[i]._mass;
  }

  const int ndiSize = mc._nuclearData->_isotopes.size();
  hipCalloc(ndiSize,isotopes);
  reactionSize = mc._nuclearData->_isotopes[0]._species[0]._reactions.size();
  const int rSizeP1 = reactionSize+1;
  assert(groupSize == mc._nuclearData->_isotopes[0]._species[0]._reactions[0]._crossSection.size());
  for (const auto &isotope : mc._nuclearData->_isotopes) {
    for (const auto &species : isotope._species) {
      assert(rSizeP1 == species._reactions.size()+1);
      for (const auto &reaction: species._reactions) {
        assert(groupSize == reaction._crossSection.size());
      }
    }
  }

  DeviceReaction *rs = nullptr;
  hipCalloc(ndiSize*rSizeP1,rs);
  double *xs = nullptr;
  hipCalloc(ndiSize*rSizeP1*groupSize,xs);
  for (int i = 0; i < ndiSize; i++) {
    isotopes[i].reactions = rs;
    for (int j = 0; j < rSizeP1; j++) {
      isotopes[i].reactions[j].crossSections = xs;
      xs += groupSize;
    }
    rs += rSizeP1;
  }

  nuBar = mc._nuclearData->_isotopes[0]._species[0]._reactions[0]._nuBar;
  for (int i = 0; i < ndiSize; i++) {
    for (int j = 1; j < rSizeP1; j++) {
      isotopes[i].reactions[j].type = mc._nuclearData->_isotopes[i]._species[0]._reactions[j-1]._reactionType;
      assert(nuBar == mc._nuclearData->_isotopes[i]._species[0]._reactions[j-1]._nuBar);
    }
    for (int k = 0; k < groupSize; k++) {
      double sum = 0;
      for (int j = 1; j < rSizeP1; j++) {
        const double xs = mc._nuclearData->_isotopes[i]._species[0]._reactions[j-1]._crossSection[k];
        sum += xs;
        isotopes[i].reactions[j].crossSections[k] = xs;
      }
      isotopes[i].reactions[0].crossSections[k] = sum;
    }
  }

  hipCalloc(PARTICLE_SIZES_SIZE,particleSizes);
  hipCalloc(TALLIES_SIZE,tallies);

  {
    const size_t vaultSize = mc._particleVaultContainer->getVaultSize();
    hipCalloc(vaultSize,processing);
    hipCalloc(vaultSize,processed);
    hipCalloc(vaultSize,extras);
  }

  {
    logLow = log(mc._nuclearData->_energies[0]);
    numGroups = mc._nuclearData->_numEnergyGroups;
    const double delta = (log(mc._nuclearData->_energies[numGroups])-logLow)/double(numGroups);
    divDelta = 1.0/delta;
  }

  timeStep = mc.time_info->time_step;
}

void Device::cycleInit(MonteCarlo &mc)
{
  const int groupSize = mc._nuclearData->_numEnergyGroups;
  const int domainSize = mc.domain.size();
  const int neighborSize = mc.domain[0].mesh._nbrRank.size();
  int cellSizeSum = 0;
  for (int i = 0; i < domainSize; i++) {
    cellSizeSum += mc.domain[i].cell_state.size();
    for (int j = 0; j < neighborSize; j++) {
      const int rank = mc.domain[i].mesh._nbrRank[j];
      domains[i].neighbors[j] = mc.particle_buffer->processor_buffer_map.at(rank);
    }
  }
  const int bytes = cellSizeSum*groupSize*sizeof(double);
  CHECK(hipMemset(domains->cells->totals,0,bytes));
  CHECK(hipMemset(domains->cells->groupTallies,0,bytes));
  CHECK(hipMemset(tallies,0,TALLIES_SIZE*sizeof(*tallies)));

  CHECK(hipMemset(particleSizes,0,PARTICLE_SIZES_SIZE*sizeof(*particleSizes)));
  const ParticleVault &vault = *(mc._particleVaultContainer->getTaskProcessingVault());
  particleSizes[PROCESSING] = vault.size();
  assert(vault.size());
  for (int i = 0; i < vault.size(); i++) processing[i] = vault[i];
}
  
void Device::cycleFinalize(MonteCarlo &mc)
{
  const int groupSize = mc._nuclearData->_numEnergyGroups;
  const int domainSize = mc.domain.size();
  for (int i = 0; i < domainSize; i++) {
    const int cellSize = mc.domain[i].cell_state.size();
    for (int j = 0; j < cellSize; j++) {
      for (int k = 0; k < groupSize; k++) {
        mc._tallies->_scalarFluxDomain[i]._task[0]._cell[j]._group[k] = domains[i].cells[j].groupTallies[k];
      }
    }
  }
  mc._tallies->_balanceTask[0]._numSegments = tallies[Tallies::SEGMENTS];
  mc._tallies->_balanceTask[0]._census = tallies[Tallies::CENSUS];
  mc._tallies->_balanceTask[0]._collision = tallies[Tallies::COLLISION];
  mc._tallies->_balanceTask[0]._escape = tallies[Tallies::ESCAPE];
  mc._tallies->_balanceTask[0]._scatter = tallies[Tallies::SCATTER];
  mc._tallies->_balanceTask[0]._absorb = tallies[Tallies::ABSORB];
  mc._tallies->_balanceTask[0]._fission = tallies[Tallies::FISSION];
  mc._tallies->_balanceTask[0]._produce = tallies[Tallies::PRODUCE];
  
  ParticleVault &processedVault = *(mc._particleVaultContainer->getTaskProcessedVault());
  processedVault.clear();
  const int processedSize = particleSizes[PROCESSED];
  MC_Base_Particle particle;
  for (int i = 0; i < processedSize; i++) {
    processed[i].set(particle);
    processedVault.pushBaseParticle(particle);
  }

  mc._particleVaultContainer->getTaskProcessingVault()->clear();
}

DeviceParticle &DeviceParticle::operator=(const MC_Base_Particle &that)
{
  coordinate.x = that.coordinate.x;
  coordinate.y = that.coordinate.y;
  coordinate.z = that.coordinate.z;
  velocity.x = that.velocity.x;
  velocity.y = that.velocity.y;
  velocity.z = that.velocity.z;
  kineticEnergy = that.kinetic_energy;
  weight = that.weight;
  timeToCensus = that.time_to_census;
  age = that.age;
  numMeanFreePaths = that.num_mean_free_paths;
  numSegments = that.num_segments;
  randomNumberSeed = that.random_number_seed;
  identifier = that.identifier;
  lastEvent = that.last_event;
  numCollisions = that.num_collisions;
  breed = that.breed;
  domain = that.domain;
  cell = that.cell;
  return *this;
}

DeviceParticle &DeviceParticle::operator=(const MC_Particle &that)
{
  return *this = MC_Base_Particle(that);
}

bool DeviceParticle::operator==(const MC_Base_Particle &that) const
{
  return ((coordinate.x == that.coordinate.x) &&
      (coordinate.y == that.coordinate.y) &&
      (coordinate.z == that.coordinate.z) &&
      (velocity.x == that.velocity.x) &&
      (velocity.y == that.velocity.y) &&
      (velocity.z == that.velocity.z) &&
      (kineticEnergy == that.kinetic_energy) &&
      (weight == that.weight) &&
      (timeToCensus == that.time_to_census) &&
      (age == that.age) &&
      (numMeanFreePaths == that.num_mean_free_paths) &&
      (numSegments == that.num_segments) &&
      (randomNumberSeed == that.random_number_seed) &&
      (identifier == that.identifier) &&
      (lastEvent == that.last_event) &&
      (numCollisions == that.num_collisions) &&
      (breed == that.breed) &&
      (domain == that.domain) &&
      (cell == that.cell));
}

void DeviceParticle::set(MC_Base_Particle &that) const
{
  that.coordinate.x = coordinate.x;
  that.coordinate.y = coordinate.y;
  that.coordinate.z = coordinate.z;
  that.velocity.x = velocity.x;
  that.velocity.y = velocity.y;
  that.velocity.z = velocity.z;
  that.kinetic_energy = kineticEnergy;
  that.weight = weight;
  that.time_to_census = timeToCensus;
  that.age = age;
  that.num_mean_free_paths = numMeanFreePaths;
  that.num_segments = numSegments;
  that.random_number_seed = randomNumberSeed;
  that.identifier = identifier;
  that.last_event = lastEvent;
  that.num_collisions = numCollisions;
  that.breed = breed;
  that.domain = domain;
  that.cell = cell;
}


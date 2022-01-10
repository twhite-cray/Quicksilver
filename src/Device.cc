#include "Device.hh"

#include "cudaUtils.hh"
#include "MonteCarlo.hh"
#include "NuclearData.hh"

void Device::init(MonteCarlo &mc)
{
  assert(domains == nullptr);
  const int domainSize = mc.domain.size();
  CHECK(hipHostMalloc(&domains,domainSize*sizeof(DeviceDomain)));

  const int gSize = mc._nuclearData->_numEnergyGroups;

  int csSizeSum = 0;
  for (int i = 0; i < domainSize; i++) csSizeSum += mc.domain[i].cell_state.size();
  
  DeviceCellState *css = nullptr;
  CHECK(hipHostMalloc(&css,csSizeSum*sizeof(DeviceCellState)));
  double *totals = nullptr;
  CHECK(hipHostMalloc(&totals,csSizeSum*gSize*sizeof(double)));
  for (int i = 0; i < domainSize; i++) {
    domains[i].cellStates = css;
    const int csSize = mc.domain[i].cell_state.size();
    for (int j = 0; j < csSize; j++) {
      css[j] = mc.domain[i].cell_state[j];
      css[j].totals = totals;
      totals += gSize;
    }
    css += csSize;
  }

  assert(mats == nullptr);
  const int matSize = mc._materialDatabase->_mat.size();
  CHECK(hipHostMalloc(&mats,matSize*sizeof(DeviceMaterial)));
  
  int isoSizeSum = 0;
  for (int i = 0; i < matSize; i++) isoSizeSum += mc._materialDatabase->_mat[i]._iso.size();

  DeviceIsotope *isos = nullptr;
  CHECK(hipHostMalloc(&isos,isoSizeSum*sizeof(DeviceIsotope)));
  for (int i = 0; i < matSize; i++) {
    mats[i].isos = isos;
    const int isoSize = mc._materialDatabase->_mat[i]._iso.size();
    mats[i].isoSize = isoSize;
    for (int j = 0; j < isoSize; j++) isos[j] = mc._materialDatabase->_mat[i]._iso[j];
    isos += isoSize;
  }

  const int ndiSize = mc._nuclearData->_isotopes.size();
  CHECK(hipHostMalloc(&isotopes,ndiSize*sizeof(DeviceNuclearDataIsotope)));
  const int rSize = mc._nuclearData->_isotopes[0]._species[0]._reactions.size()+1;
  assert(gSize == mc._nuclearData->_isotopes[0]._species[0]._reactions[0]._crossSection.size());
  for (const auto &isotope : mc._nuclearData->_isotopes) {
    for (const auto &species : isotope._species) {
      assert(rSize == species._reactions.size()+1);
      for (const auto &reaction: species._reactions) {
        assert(gSize == reaction._crossSection.size());
      }
    }
  }

  DeviceReaction *rs = nullptr;
  CHECK(hipHostMalloc(&rs,ndiSize*rSize*sizeof(DeviceReaction)));
  double *xs = nullptr;
  CHECK(hipHostMalloc(&xs,ndiSize*rSize*gSize*sizeof(double)));
  for (int i = 0; i < ndiSize; i++) {
    isotopes[i].reactions = rs;
    for (int j = 0; j < rSize; j++) {
      isotopes[i].reactions[j].crossSections = xs;
      xs += gSize;
    }
    rs += rSize;
  }

  for (int i = 0; i < ndiSize; i++) {
    for (int k = 0; k < gSize; k++) {
      double sum = 0;
      for (int j = 1; j < rSize; j++) {
        const double xs = mc._nuclearData->_isotopes[i]._species[0]._reactions[j-1]._crossSection[k];
        sum += xs;
        isotopes[i].reactions[j].crossSections[k] = xs;
      }
      isotopes[i].reactions[0].crossSections[k] = sum;
    }
  }
}

void Device::cycleInit(MonteCarlo &mc)
{
  const int gSize = mc._nuclearData->_numEnergyGroups;
  const int domainSize = mc.domain.size();
  int csSizeSum = 0;
  for (int i = 0; i < domainSize; i++) csSizeSum += mc.domain[i].cell_state.size();
  memset(domains->cellStates->totals,0,csSizeSum*gSize*sizeof(double));
}
  


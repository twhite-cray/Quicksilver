#include "Device.hh"

#include "cudaUtils.hh"
#include "MonteCarlo.hh"
#include "NuclearData.hh"

void Device::init(MonteCarlo &mc)
{
  assert(domain == nullptr);
  const int domainSize = mc.domain.size();
  CHECK(hipHostMalloc(&domain,domainSize*sizeof(DeviceDomain)));

  int csSizeSum = 0;
  for (int i = 0; i < domainSize; i++) csSizeSum += mc.domain[i].cell_state.size();
  
  DeviceCellState *cs = nullptr;
  CHECK(hipHostMalloc(&cs,csSizeSum*sizeof(DeviceCellState)));
  for (int i = 0; i < domainSize; i++) {
    domain[i].cellState = cs;
    const int csSize = mc.domain[i].cell_state.size();
    for (int j = 0; j < csSize; j++) cs[j] = mc.domain[i].cell_state[j];
    cs += csSize;
  }

  assert(mat == nullptr);
  const int matSize = mc._materialDatabase->_mat.size();
  CHECK(hipHostMalloc(&mat,matSize*sizeof(DeviceMaterial)));
  
  int isoSizeSum = 0;
  for (int i = 0; i < matSize; i++) isoSizeSum += mc._materialDatabase->_mat[i]._iso.size();

  DeviceIsotope *iso = nullptr;
  CHECK(hipHostMalloc(&iso,isoSizeSum*sizeof(DeviceIsotope)));
  for (int i = 0; i < matSize; i++) {
    mat[i].iso = iso;
    const int isoSize = mc._materialDatabase->_mat[i]._iso.size();
    for (int j = 0; j < isoSize; j++) iso[j] = mc._materialDatabase->_mat[i]._iso[j];
    iso += isoSize;
  }

  const int ndiSize = mc._nuclearData->_isotopes.size();
  CHECK(hipHostMalloc(&isotopes,ndiSize*sizeof(DeviceNuclearDataIsotope)));
  const int groupSize = mc._nuclearData->_isotopes[0]._species[0]._reactions[0]._crossSection.size();
  for (const auto &isotope : mc._nuclearData->_isotopes) {
    for (const auto &species : isotope._species) {
      for (const auto &reaction: species._reactions) {
        assert(groupSize == reaction._crossSection.size());
      }
    }
  }

  double *tcs = nullptr;
  CHECK(hipHostMalloc(&tcs,ndiSize*groupSize*sizeof(double)));
  for (int i = 0; i < ndiSize; i++) {
    isotopes[i].totalCrossSections = tcs;
    const auto &reactions = mc._nuclearData->_isotopes[i]._species[0]._reactions;
    for (int j = 0; j < groupSize; j++) {
      double sum = 0;
      for (const auto &reaction : reactions) sum += reaction._crossSection[j];
      tcs[j] = sum;
    }
    tcs += groupSize;
  }
}


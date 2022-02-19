#ifndef MACROSCOPIC_CROSS_SECTION_HH
#define MACROSCOPIC_CROSS_SECTION_HH

#include "DeclareMacro.hh"
#include "Device.hh"

//----------------------------------------------------------------------------------------------------------------------
//  Routine MacroscopicCrossSection calculates the number-density-weighted macroscopic cross
//  section of a cell.
//
//  A reactionIndex of -1 means total cross section.
//----------------------------------------------------------------------------------------------------------------------
__host__ __device__ static inline double macroscopicCrossSection(const Device &device, int reactionIndex, int domainIndex, int cellIndex,
                               int isoIndex, int energyGroup)
{
   // The cell number density is the fraction of the atoms in cell
   // volume of this isotope.  We set this (elsewhere) to 1/nIsotopes.
   // This is a statement that we treat materials as if all of their
   // isotopes are present in equal amounts

   const DeviceCell &cell = device.domains[domainIndex].cells[cellIndex];
   const double cellNumberDensity = cell.cellNumberDensity;
   const int globalMatIndex = cell.material;

   const DeviceIsotope &iso = device.mats[globalMatIndex].isos[isoIndex];
   const double atomFraction = iso.atomFraction;
   const int isotopeGid = iso.gid;

   if ( atomFraction == 0.0 || cellNumberDensity == 0.0) { return 1e-20; }

   const int rix = (reactionIndex < 0) ? 0 : reactionIndex+1;
   // 0 -> total cross section
   const double microscopicCrossSection = device.isotopes[isotopeGid].reactions[rix].crossSections[energyGroup];

   return atomFraction * cellNumberDensity * microscopicCrossSection;
}


//----------------------------------------------------------------------------------------------------------------------
//  Routine weightedMacroscopicCrossSection calculates the number-density-weighted
//  macroscopic cross section of the collection of isotopes in a cell.
//dfr Weighted is a bit of a misnomer here, since there is no weighting
//applied by this routine.  In Mercury we would weight for multiple
//materials in a cell.
//----------------------------------------------------------------------------------------------------------------------
__host__ __device__ static inline double weightedMacroscopicCrossSection(Device &device, int taskIndex, int domainIndex,
                                       int cellIndex, int energyGroup)
{
   const double xs = device.domains[domainIndex].cells[cellIndex].totals[energyGroup];
   if (xs > 0.0) return xs;

   const int globalMatIndex = device.domains[domainIndex].cells[cellIndex].material;
   const int nIsotopes = device.mats[globalMatIndex].isoSize;
   double sum = 0.0;
   for (int isoIndex = 0; isoIndex < nIsotopes; isoIndex++)
   {
      sum += macroscopicCrossSection(device, -1, domainIndex, cellIndex, isoIndex, energyGroup);
   }
   device.domains[domainIndex].cells[cellIndex].totals[energyGroup] = sum;
   return sum;
}

#endif

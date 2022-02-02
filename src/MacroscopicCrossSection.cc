#include "MacroscopicCrossSection.hh"
#include "MonteCarlo.hh"
#include "MaterialDatabase.hh"
#include "NuclearData.hh"
#include "MC_Cell_State.hh"
#include "DeclareMacro.hh"

//----------------------------------------------------------------------------------------------------------------------
//  Routine MacroscopicCrossSection calculates the number-density-weighted macroscopic cross
//  section of a cell.
//
//  A reactionIndex of -1 means total cross section.
//----------------------------------------------------------------------------------------------------------------------
double macroscopicCrossSection(const Device &device, int reactionIndex, int domainIndex, int cellIndex,
                               int isoIndex, int energyGroup)
{
   // The cell number density is the fraction of the atoms in cell
   // volume of this isotope.  We set this (elsewhere) to 1/nIsotopes.
   // This is a statement that we treat materials as if all of their
   // isotopes are present in equal amounts

   const DeviceCellState &cellState = device.domains[domainIndex].cellStates[cellIndex];
   const double cellNumberDensity = cellState.cellNumberDensity;
   const int globalMatIndex = cellState.material;

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
double weightedMacroscopicCrossSection(Device &device, int taskIndex, int domainIndex,
                                       int cellIndex, int energyGroup)
{
   const double xs = device.domains[domainIndex].cellStates[cellIndex].totals[energyGroup];
   if (xs > 0.0) return xs;

   const int globalMatIndex = device.domains[domainIndex].cellStates[cellIndex].material;
   const int nIsotopes = device.mats[globalMatIndex].isoSize;
   double sum = 0.0;
   for (int isoIndex = 0; isoIndex < nIsotopes; isoIndex++)
   {
      sum += macroscopicCrossSection(device, -1, domainIndex, cellIndex, isoIndex, energyGroup);
   }
   device.domains[domainIndex].cellStates[cellIndex].totals[energyGroup] = sum;
   return sum;
}

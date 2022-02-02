#ifndef MACROSCOPIC_CROSS_SECTION_HH
#define MACROSCOPIC_CROSS_SECTION_HH

#include "DeclareMacro.hh"

struct Device;

double macroscopicCrossSection(const Device &device, int reactionIndex, int domainIndex, int cellIndex,
                               int isoIndex, int energyGroup);

double weightedMacroscopicCrossSection(Device &device, int taskIndex, int domainIndex,
                                       int cellIndex, int energyGroup);

#endif

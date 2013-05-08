#ifndef GENERATE_COMBOS_H
#define GENERATE_COMBOS_H

#include "Setup.h"

//This function generates combinations after dependency checks have identifiied independent
//input output pairs. The process of generation includes generating binary bit vectors and
//metbaolite coefficeints of the combinations.
//Param metabolite: The metabolite being balanced
//Param numberOfBin: Number of bins in combination bin array
//Param nextFreePathwayIndex: Index in the pathway buffers to store new pathways
int generateCombinations(int metabolite, int numberOfBins, int nextFreePathwayIndex);

#endif

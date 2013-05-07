#ifndef SORT_PATHWAYS_H
#define SORT_PATHWAYS_H

//TODO: define CPU code the calls kernel
#include "Setup.h"

void sortInputsOutputs(float *d_metaboliteCoefficients, int pathwayCounts,
        BinaryVector *d_reactions, int metaboliteCount, int numInputs,
        int numOutputs, int metaboliteToRemove);
void dependencyCheck(int numInputs, int numOutputs, int batch_number);

#endif

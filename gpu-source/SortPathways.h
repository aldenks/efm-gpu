#ifndef SORT_PATHWAYS_H
#define SORT_PATHWAYS_H

//TODO: define CPU code the calls kernel
#include "Setup.h"

#define MAX_THREADS_PER_BLOCK 1024

void sortInputsOutputs(float *d_metaboliteCoefficients, int pathwayCounts, 
		       BinaryVector *d_reactions, int metaboliteCount, int numInputs, 
		       int numOutputs, int metaboliteToRemove);

#endif

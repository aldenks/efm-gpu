#ifndef SETUP_H
#define	SETUP_H

#include "Network.h"

//Tolerence for positive zero
#define ZERO 1e-10
//Tolerence for negative zero
#define NEG_ZERO -ZERO

//Maximum number of pathways that can be stored in the memory
#define MAX_PATHWAYS 1000000

// number of combinations tested per batch
extern int batchSize;

//Binary Vector Datatype
typedef unsigned long BinaryVector;

//Metabolic network
extern Network network;

//Device data for binary vector
extern BinaryVector* d_binaryVectors;
//Device data for metabolite coefficients
//2D array. Each row represents to the coefficients for each metabolite in a pathway
extern float* d_metaboliteCoefficients;
//Device data for combinations
extern int* d_combinationBins;
//Device flags for balanced metabolites. A true means the metabolite is balanced
extern bool* d_balancedMetabolites;
//Device data for count of input pathways for each metabolite
extern int* d_metaboliteInputPathwayCounts;
//Device data for count of output pathways for each metabolite
extern int* d_metaboliteOutputPathwayCounts;

//Host data for count of input pathways for each metabolite
extern int* h_metaboliteInputPathwayCounts;
//Host data for count of output pathways for each metabolite
extern int* h_metaboliteOutputPathwayCounts;

//Number of remaining metabolites to be balanced
extern int remainingMetabolites;
//Number of metabolites
extern int metaboliteCount;
//Number of current pathways
extern int pathwayCount;

// Bins for each thread's newly found independent pathways
//  bins are organized in column major for memory coalescing
#define BIN_MAX_ENTRIES 13107200 // 50mb of 4 byte ints
extern int *d_newPathwayBins;
extern int *d_newPathwayBinCounts;
extern int *h_newPathwayBinCounts;
// Indices into the binaryVectors & metaboliteCoefficients for
//  threads to begin writing out to
extern int *d_newPathwayWriteIndices;
extern int *h_newPathwayWriteIndices;

//Initializes the network and gpu memory
//Ret: returns true if memory was succesfully allocated, false otherwise
bool setup();
//Frees the allocated memory
void freeResources();


#endif	/* SETUP_H */
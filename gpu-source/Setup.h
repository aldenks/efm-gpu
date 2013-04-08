#ifndef SETUP_H
#define	SETUP_H

#include "Network.h"

//Tolerence for positive zero
#define ZERO 1e-10
//Tolerence for negative zero
#define NEG_ZERO -ZERO

//Maximum number of pathways that can be stored in the memory
#define MAX_PATHWAYS 1000000

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
extern bool* d_balancedMetabolies;

//Host data for binary vector
extern BinaryVector* h_binaryVectors;
//Number of remaining metabolites to be balanced
extern int remainingMetabolites;
//Number of metbaolites
extern int metaboliteCount;
//Number of current pathways
extern int pathwayCount;

//Initializes the network and gpu memory
//Ret: returns true if memory was succesfully allocated, false otherwise
bool setup();
//Frees the allocated memory
void freeResources();


#endif	/* SETUP_H */

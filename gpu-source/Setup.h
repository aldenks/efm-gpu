#ifndef SETUP_H
#define	SETUP_H

#include "Network.h"

//Tolerence for positive zero
#define ZERO 10e-10
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
extern float* d_metaboliteCoefficients;
//Device data for combinations
extern int* d_combinationBins;

//Host data for binary vector
extern BinaryVector* h_binaryVectors;

//Initializes the network and gpu memory
//Ret: returns true if memory was succesfully allocated, false otherwise
bool setup();
//Frees the allocated memory
void freeResources();


#endif	/* SETUP_H */

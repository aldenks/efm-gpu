#include "Setup.h"

//Metabolic network
Network network;

//Device data for binary vector
BinaryVector* d_binaryVectors;
//Device data for metabolite coefficients
//2D array. Each row represents to the coefficients for each metabolite in a pathway
float* d_metaboliteCoefficients;
//Device data for combinations
int* d_combinationBins;
//Device flags for balanced metabolites. A true means the metabolite is balanced
bool* d_balancedMetabolites;

//Host data for binary vector
BinaryVector* h_binaryVectors;
//Number of remaining metabolites to be balanced
int remainingMetabolites;
//Number of metbaolites
int metaboliteCount;
//Number of pathways
int pathwayCount;

bool allocateMemory() {
   //Allocating memory for host binary vectors
   h_binaryVectors = (BinaryVector*) malloc(MAX_PATHWAYS * sizeof (BinaryVector));
   if (!h_binaryVectors) {
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, "Unable to allocate memory for h_binaryVectors");
      return false;
   }
   cudaError error;

   //Allocating memory for device binary vectors
   d_binaryVectors = NULL;
   error = cudaMalloc((void**) &d_binaryVectors, MAX_PATHWAYS * sizeof (BinaryVector));
   if (error != cudaSuccess) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for d_binaryVectors\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));
      return false;
   }

   //Allocating memory for device metabolite coefficients
   d_metaboliteCoefficients = NULL;
   error = cudaMalloc((void**) &d_metaboliteCoefficients, MAX_PATHWAYS * network.metabolites.size() * sizeof (float));
   if (error != cudaSuccess) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for d_metaboliteCoefficients\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));
      return false;
   }

   //Allocating memory for device metabolite balanced state
   d_balancedMetabolites = NULL;
   error = cudaMalloc((void**) &d_balancedMetabolites, network.metabolites.size() * sizeof(bool));
   if (error != cudaSuccess) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for d_balancedMetabolites\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));
      return false;
   }

   return true;
}

//Generates bit packed representation of reactions

BinaryVector packReaction(int reaction) {
   return 1 << reaction;
}

//Initializes the network and gpu memory
//Ret: returns true if memory was succesfully allocated, false otherwise

bool setup() {
   //Allocate memory
   if (!allocateMemory()) {
      return false;
   }
   int metCoeffSize = network.reactions.size() * network.metabolites.size() * sizeof (float);
   float* metCoeff = (float*) malloc(metCoeffSize);
   //Initialize cpu memory
   for (int r = 0, mIndex = 0; r < network.reactions.size(); r++) {
      h_binaryVectors[r] = packReaction(r);
      for (int m = 0; m < network.metabolites.size(); m++) {
         metCoeff[mIndex++] = network.s[m][r];
      }
   }
   //Initialize gpu memory
   cudaMemcpy(d_binaryVectors, h_binaryVectors, network.reactions.size() * sizeof (BinaryVector), cudaMemcpyHostToDevice);
   cudaMemcpy(d_metaboliteCoefficients, metCoeff, metCoeffSize, cudaMemcpyHostToDevice);
   free(metCoeff);

   // initialize counts
   // TODO Ehsan, are these right?
   pathwayCount = network.reactions.size();
   metaboliteCount = network.metabolites.size();
   remainingMetabolites = metaboliteCount;

   // intialize all metabolites balanced state to false
   cudaMemset(d_balancedMetabolites, 0, metaboliteCount * sizeof(bool));

   return true;
}

//Frees the allocated memory

void freeResources() {
   if (h_binaryVectors) {
      free(h_binaryVectors);
   }
   if (d_binaryVectors) {
      cudaFree(d_binaryVectors);
   }
   if (d_metaboliteCoefficients) {
      cudaFree(d_metaboliteCoefficients);
   }
   if (d_combinationBins) {
      cudaFree(d_combinationBins);
   }
}

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
//Device data for count of input pathways for each metabolite
int* d_metaboliteInputPathwayCounts;
//Device data for count of output pathways for each metabolite
int* d_metaboliteOutputPathwayCounts;

//Host data for count of input pathways for each metabolite
int* h_metaboliteInputPathwayCounts;
//Host data for count of output pathways for each metabolite
int* h_metaboliteOutputPathwayCounts;

//Number of remaining metabolites to be balanced
int remainingMetabolites;
//Number of metabolites
int metaboliteCount;
//Number of current pathways
int pathwayCount;

bool allocateMemory() {
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
   error = cudaMalloc((void**) &d_metaboliteCoefficients, MAX_PATHWAYS * metaboliteCount * sizeof (float));
   if (error != cudaSuccess) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for d_metaboliteCoefficients\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));
      return false;
   }

   //Allocating memory for device metabolite balanced state
   d_balancedMetabolites = NULL;
   error = cudaMalloc((void**) &d_balancedMetabolites, metaboliteCount * sizeof (bool));
   if (error != cudaSuccess) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for d_balancedMetabolites\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));
      return false;
   }

   //Device data for count of input pathways for each metabolite
   d_metaboliteInputPathwayCounts = NULL;
   error = cudaMalloc((void**) &d_metaboliteInputPathwayCounts, metaboliteCount * sizeof (int));
   if (error != cudaSuccess) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for d_metaboliteInputPathwayCounts\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));
      return false;
   }

   //Device data for count of input pathways for each metabolite
   d_metaboliteOutputPathwayCounts = NULL;
   error = cudaMalloc((void**) &d_metaboliteOutputPathwayCounts, metaboliteCount * sizeof (int));
   if (error != cudaSuccess) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for d_metaboliteOutputPathwayCounts\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));
      return false;
   }

   //Host data for count of input pathways for each metabolite
   h_metaboliteInputPathwayCounts = (int*) malloc(metaboliteCount * sizeof (int));
   if (!h_metaboliteInputPathwayCounts) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for h_metaboliteInputPathwayCounts\n");
      return false;
   }

   //Host data for count of output pathways for each metabolite
   h_metaboliteOutputPathwayCounts = (int*) malloc(metaboliteCount * sizeof (int));
   if (!h_metaboliteInputPathwayCounts) {
      fprintf(stderr, "Setup.cu:allocateMemory() Unable to allocate memory for h_metaboliteOutputPathwayCounts\n");
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
   //Initialize pathway count. It is equal to the number of reaction in the network
   //after splitting the reversible reactions
   pathwayCount = network.reactions.size();

   //Initialize metabolite count. This is equal to the number of internal metabolites
   metaboliteCount = 0;
   for (int i = 0; i < network.metabolites.size(); i++) {
      if (!network.external[i]) {
         metaboliteCount++;
      }
   }
   //The number of remaining metabolites is equal to the number of internal metabolites
   //at the beginning of the algorithm
   remainingMetabolites = metaboliteCount;

   //Allocate memory
   if (!allocateMemory()) {
      return false;
   }

   //Allocating memory for host binary vectors
   int binaryVectorSize = pathwayCount * sizeof (BinaryVector);
   BinaryVector* h_binaryVectors = (BinaryVector*) malloc(binaryVectorSize);
   int metCoeffSize = pathwayCount * metaboliteCount * sizeof (float);
   float* metCoeff = (float*) malloc(metCoeffSize);
   //Initialize cpu memory
   for (int r = 0, mIndex = 0; r < pathwayCount; r++) {
      h_binaryVectors[r] = packReaction(r);
      for (int m = 0; m < network.metabolites.size(); m++) {
         if (!network.external[m]) {
            metCoeff[mIndex++] = network.s[m][r];
         }
      }
   }
   //Initialize gpu memory
   cudaMemcpy(d_binaryVectors, h_binaryVectors, binaryVectorSize, cudaMemcpyHostToDevice);
   cudaMemcpy(d_metaboliteCoefficients, metCoeff, metCoeffSize, cudaMemcpyHostToDevice);

   //Free cpu memory used for initialization of gpu memory
   free(h_binaryVectors);
   free(metCoeff);

   // intialize all metabolites balanced state to false
   cudaMemset(d_balancedMetabolites, 0, metaboliteCount * sizeof (bool));

   return true;
}

//Frees the allocated memory

void freeResources() {
   if (d_binaryVectors) {
      cudaFree(d_binaryVectors);
   }
   if (d_metaboliteCoefficients) {
      cudaFree(d_metaboliteCoefficients);
   }
   if (d_combinationBins) {
      cudaFree(d_combinationBins);
   }
   if (d_balancedMetabolites) {
      cudaFree(d_balancedMetabolites);
   }
   if (d_metaboliteInputPathwayCounts) {
      cudaFree(d_metaboliteInputPathwayCounts);
   }
   if (d_metaboliteOutputPathwayCounts) {
      cudaFree(d_metaboliteOutputPathwayCounts);
   }
   if (h_metaboliteInputPathwayCounts) {
      free(h_metaboliteInputPathwayCounts);
   }
   if (h_metaboliteOutputPathwayCounts) {
      free(h_metaboliteOutputPathwayCounts);
   }
}

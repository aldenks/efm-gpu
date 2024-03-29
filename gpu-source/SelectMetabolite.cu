#include "SelectMetabolite.h"
#include "Setup.h"

#define BALANCED_MARKER -1

// Kernel with one thread per metabolite. Sums input/output counts and writes
// them to the metabolite's index in inputCounts and outputCounts.
__global__
void computeMetaboliteInputOutputCounts(float* metaboliteCoefficients, int pathwayStartIndex, int pathwayCount, int metaboliteCount, bool* balancedMetabolites, int* inputCounts, int* outputCounts) {

   int inputCount = 0;
   int outputCount = 0;
   int m = blockIdx.x * blockDim.x + threadIdx.x;
   if (m >= metaboliteCount) {
      return;
   }

   int end_i = pathwayStartIndex + pathwayCount;
   for (int p = pathwayStartIndex; p < end_i; p++) {
      float coeff = metaboliteCoefficients[ circularIndex(p * metaboliteCount + m) ];
      if (coeff > ZERO) {
         outputCount++;
      } else if (coeff < NEG_ZERO) {
         inputCount++;
      }
   }

   inputCounts[m] = inputCount;
   outputCounts[m] = outputCount;
   // insert sentinel value if the metabolite is already balanced
   if (balancedMetabolites[m]) {
      inputCounts[m] = BALANCED_MARKER;
      outputCounts[m] = BALANCED_MARKER;
   }
}

// Returns the index of the next metabolite to remove.
// Side effect: h_{input|output}Counts contain the correct counts
//  for each metabolite.
int getNextMetabolite(float* d_metaboliteCoefficients, int pathwayStartIndex, int pathwayCount, int metaboliteCount, bool* d_balancedMetabolites, int* d_inputCounts, int* d_outputCounts, int* h_inputCounts, int* h_outputCounts) {

   int threads_per_count_block = 256;
   int blocks_per_grid = ceil(((float) metaboliteCount) / threads_per_count_block);
   computeMetaboliteInputOutputCounts << <blocks_per_grid, threads_per_count_block >> >
           (d_metaboliteCoefficients, pathwayStartIndex, pathwayCount, metaboliteCount, d_balancedMetabolites, d_inputCounts, d_outputCounts);
   int count_mem_size = metaboliteCount * sizeof (int);
   cudaMemcpy(h_inputCounts, d_inputCounts, count_mem_size, cudaMemcpyDeviceToHost);
   cudaMemcpy(h_outputCounts, d_outputCounts, count_mem_size, cudaMemcpyDeviceToHost);

   // select the metabolite with the minimum product of inputs and outputs
   int min_i = 0;
   long min_product = LONG_MAX;
   for (int i = 0; i < metaboliteCount; i++) {
      if (h_inputCounts[i] == BALANCED_MARKER) {
         continue;
      }
      int product = h_inputCounts[i] * h_outputCounts[i];
      if (product < min_product) {
         min_i = i;
         min_product = product;
         if (min_product == 0) return min_i; // stop if we've already found a "best"
      }
   }
   return min_i;
}

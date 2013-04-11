#include "SelectMetabolite.h"

#define BALANCED_MARKER -1

__global__
void computeMetaboliteInputOutputCounts(float* metaboliteCoefficients, int pathwayCount, int metaboliteCount, bool* balancedMetabolites, int* inputCounts, int* outputCounts) {

   int inputCount = 0;
   int outputCount = 0;
   int m = blockIdx.x * blockDim.x + threadIdx.x;
   if (m >= metaboliteCount) {
      return;
   }

   for (int p = 0; p < pathwayCount; p++) {
      float coeff = metaboliteCoefficients[ p * metaboliteCount + m ];
      if (coeff > 0) {
         outputCount++;
      } else if (coeff < 0) {
         inputCount++;
      }
   }

   inputCounts[m] = inputCount;
   outputCounts[m] = outputCount;
   // insert sentinel value if the metabolite is already balanced
   if (balancedMetabolites[m]) {
      inputCounts[m]  = BALANCED_MARKER;
      outputCounts[m] = BALANCED_MARKER;
   }
}

int getNextMetabolite(float* d_metaboliteCoefficients, int pathwayCount, int metaboliteCount, bool* d_balancedMetabolites, int* d_inputCounts, int* d_outputCounts, int* h_inputCounts, int* h_outputCounts) {

   int threads_per_count_block = 256;
   int blocks_per_grid = ceil(((double) metaboliteCount) / threads_per_count_block);
   computeMetaboliteInputOutputCounts << <blocks_per_grid, threads_per_count_block >> >
           (d_metaboliteCoefficients, pathwayCount, metaboliteCount, d_balancedMetabolites,
           d_inputCounts, d_outputCounts);
   int count_mem_size = metaboliteCount * sizeof (int);
   cudaMemcpy(h_inputCounts, d_inputCounts, count_mem_size, cudaMemcpyDeviceToHost);
   cudaMemcpy(h_outputCounts, d_outputCounts, count_mem_size, cudaMemcpyDeviceToHost);

   // select the metabolite with the minimum product of inputs and outputs
   int min_i;
   long min_product = LONG_MAX;
   for (int i = 0; i < metaboliteCount; i++) {
      if (h_inputCounts[i] == BALANCED_MARKER) continue;
      int product = h_inputCounts[i] * h_inputCounts[i];
      if (product < min_product) {
         min_i = i;
         min_product = product;
         if (min_product == 0) return min_i; // stop if we've already found a "best"
      }
   }
   return min_i;
}

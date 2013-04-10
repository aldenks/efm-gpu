#include "EFMGenerator.h"
#include <iostream>
using namespace std;

//input and output counts for each metabolite. holds INT_MAX if metabolite is already processed
int* h_inputCounts;
int* h_outputCounts;
int* d_inputCounts;
int* d_outputCounts;

//Generates EFMs
void generateEFMs() {

   // allocate memory to hold input/output counts
   int count_mem_size = metaboliteCount * sizeof(int);
   h_inputCounts  = (int*)malloc(count_mem_size);
   h_outputCounts = (int*)malloc(count_mem_size);
   cudaMalloc((void**) &d_inputCounts,  count_mem_size);
   cudaMalloc((void**) &d_outputCounts, count_mem_size);

   while (remainingMetabolites > 0) {
      // find metabolite to remove
      int metabolite = getNextMetabolite(d_metaboliteCoefficients, pathwayCount, metaboliteCount,
                                         d_balancedMetabolites, d_inputCounts, d_outputCounts,
                                         h_inputCounts, h_outputCounts);

      //Copy the bit vectors from gpu to cpu
      //Sort the bit vectors for inputs outputs and non-participating pathways
      //Copy the bit vectors from cpu to gpu
      //Setup bins
      //Call kernel to generate combinations
      remainingMetabolites--;
   }
}

#include "EFMGenerator.h"

//Generates EFMs

void generateEFMs() {

   while (remainingMetabolites > 0) {
      // find metabolite to remove
      int metabolite = getNextMetabolite(d_metaboliteCoefficients, pathwayCount, metaboliteCount,
              d_balancedMetabolites, d_metaboliteInputPathwayCounts, d_metaboliteOutputPathwayCounts,
              h_metaboliteInputPathwayCounts, h_metaboliteOutputPathwayCounts);

      //Copy the bit vectors from gpu to cpu
      //Sort the bit vectors for inputs outputs and non-participating pathways
      sortInputsOutputs(d_metaboliteCoefficients, pathwayCount, d_binaryVectors,
			  metaboliteCount, h_metaboliteInputPathwayCounts[metabolite],
			  h_metaboliteOutputPathwayCounts[metabolite], metabolite);
      int divisor = h_metaboliteInputPathwayCounts[metabolite] > 0 ? 
                        h_metaboliteInputPathwayCounts[metabolite] : 1;
      batchSize = BIN_MAX_ENTRIES / divisor; // int division
      //Call kernel to generate combinations

      markMetaboliteBalanced(metabolite, d_balancedMetabolites);
      remainingMetabolites--;
   }
}

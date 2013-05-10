#include "EFMGenerator.h"
#include "SortPathways.h"
#include "GenerateCombos.h"

//Generates EFMs

void generateEFMs() {
   int serial = 1;
   printf("Metabolites=%d, Reactions=%d\n", metaboliteCount, pathwayCount);
   while (remainingMetabolites > 0) {
      // find metabolite to remove
      int metabolite = getNextMetabolite(d_metaboliteCoefficients, pathwayStartIndex,
              pathwayCount, metaboliteCount, d_balancedMetabolites,
              d_metaboliteInputPathwayCounts, d_metaboliteOutputPathwayCounts,
              h_metaboliteInputPathwayCounts, h_metaboliteOutputPathwayCounts);
      printf("%2d %2d %6s in=%-2d out=%-2d combinations=%-10d\n", serial++, metabolite, network.metabolites[metabolite].c_str(), h_metaboliteInputPathwayCounts[metabolite], h_metaboliteOutputPathwayCounts[metabolite], h_metaboliteInputPathwayCounts[metabolite] * h_metaboliteOutputPathwayCounts[metabolite]);

      //Sort the bit vectors for inputs outputs and non-participating pathways
      sortInputsOutputs(h_metaboliteInputPathwayCounts[metabolite], h_metaboliteOutputPathwayCounts[metabolite], metabolite);
      //Variable storing count for new combinations
      int newPathwayCount = 0;
      //Number of batches for dependency checking

      //Index in pathway buffers where new combinations are stored
      int nextFreePathwayIndex = circularIndex(pathwayStartIndex + pathwayCount);

      // As long as we have inputs or outputs
      if (h_metaboliteInputPathwayCounts[metabolite] != 0 && h_metaboliteOutputPathwayCounts[metabolite] != 0) {
         //Prevent division by 0
         int divisor = h_metaboliteInputPathwayCounts[metabolite] > 0 ? h_metaboliteInputPathwayCounts[metabolite] : 1;
         batchSize = BIN_MAX_ENTRIES / divisor; // int division
	 int num_batches = ((h_metaboliteOutputPathwayCounts[metabolite] - 1) / batchSize) + 1;

         printf("PathwayCount=%-5d PathwayStartIndex=%-5d NextFreeIndex=%-5d Batches=%-2d\n", pathwayCount, pathwayStartIndex, nextFreePathwayIndex, num_batches);

         //Generate combinations and check dependencies
         for (int i = 0; i < num_batches; ++i) {
            dependencyCheck(h_metaboliteInputPathwayCounts[metabolite], h_metaboliteOutputPathwayCounts[metabolite], i);
            newPathwayCount += generateCombinations(metabolite, h_metaboliteInputPathwayCounts[metabolite], circularIndex(nextFreePathwayIndex + newPathwayCount));
         }
      }

      // "Throw away" inputs and outputs by updating starting index and count of pathways
      pathwayStartIndex = circularIndex(pathwayStartIndex + h_metaboliteInputPathwayCounts[metabolite] + h_metaboliteOutputPathwayCounts[metabolite]);

      //Update the number of pathways
      pathwayCount += newPathwayCount - h_metaboliteInputPathwayCounts[metabolite] - h_metaboliteOutputPathwayCounts[metabolite];

      // Mark this metabolite as complete
      markMetaboliteBalanced(metabolite, d_balancedMetabolites);
      remainingMetabolites--;

      printf("NewPathwayCount=%-5d PathwayCount=%-5d\n", newPathwayCount, pathwayCount);
      getchar();
   }
}

#include "EFMGenerator.h"
#include "SortPathways.h"
#include "GenerateCombos.h"

//Generates EFMs

void generateEFMs() {
   int serial = 1;
   printf("Metabolites=%d, Reactions=%d\n",metaboliteCount, pathwayCount);
   while (remainingMetabolites > 0) {
      // find metabolite to remove
      int metabolite = getNextMetabolite(d_metaboliteCoefficients, pathwayStartIndex,
              pathwayCount, metaboliteCount, d_balancedMetabolites,
              d_metaboliteInputPathwayCounts, d_metaboliteOutputPathwayCounts,
              h_metaboliteInputPathwayCounts, h_metaboliteOutputPathwayCounts);
      printf("%2d %2d %6s in=%-2d out=%-2d combinations=%-10d\n", serial++, metabolite, network.metabolites[metabolite].c_str(), h_metaboliteInputPathwayCounts[metabolite], h_metaboliteOutputPathwayCounts[metabolite], h_metaboliteInputPathwayCounts[metabolite] * h_metaboliteOutputPathwayCounts[metabolite]);
      //Copy the bit vectors from gpu to cpu
      //Sort the bit vectors for inputs outputs and non-participating pathways
      sortInputsOutputs(h_metaboliteInputPathwayCounts[metabolite], h_metaboliteOutputPathwayCounts[metabolite], metabolite);
      int divisor = h_metaboliteInputPathwayCounts[metabolite] > 0 ? h_metaboliteInputPathwayCounts[metabolite] : 1;
      batchSize = BIN_MAX_ENTRIES / divisor; // int division
      //Call kernel to generate combinations
      int num_batches = (h_metaboliteOutputPathwayCounts[metabolite] / pathwayCount) + 1;
      int nextFreePathwayIndex = pathwayStartIndex + pathwayCount;
      int newPathwayCount = 0;
      printf("PathwayCount=%-5d PathwayStartIndex=%-5d NextFreeIndex=%-5d Batches=%-2d\n", pathwayCount, pathwayStartIndex, nextFreePathwayIndex, num_batches);
      for (int i = 0; i < num_batches; ++i) {
         dependencyCheck(h_metaboliteInputPathwayCounts[metabolite], h_metaboliteOutputPathwayCounts[metabolite], i);
         newPathwayCount += generateCombinations(metabolite, h_metaboliteInputPathwayCounts[metabolite], circularIndex(nextFreePathwayIndex + newPathwayCount));
      }
      pathwayStartIndex = circularIndex(pathwayStartIndex + h_metaboliteInputPathwayCounts[metabolite] + h_metaboliteOutputPathwayCounts[metabolite]);
      pathwayCount += newPathwayCount - h_metaboliteInputPathwayCounts[metabolite] - h_metaboliteOutputPathwayCounts[metabolite];
      markMetaboliteBalanced(metabolite, d_balancedMetabolites);
      remainingMetabolites--;
      printf("NewPathwayCount=%-5d\n",newPathwayCount);
      printf("PathwayCount=%-5d PathwayStartIndex=%-5d NextFreeIndex=%-5d Batches=%-2d\n", pathwayCount, pathwayStartIndex, nextFreePathwayIndex, num_batches);
      getchar();
   }
}

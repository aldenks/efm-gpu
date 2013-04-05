#include <stdio.h>
#include <stdlib.h>

#include "EFMGenerator.h"

Network network;

void execute(const char* file) {
   if (!network.readNetworkFile(file)) {
      printf("Error loading network file\n");
      return;
   }
   numReactions = network.reactions.size();
   reactionBitDataLength = ((numReactions - 1) / (8 * sizeof (int))) + 1;
   numMetabolitesRemaining = numMetabolites = network.metabolites.size();
   Pathway* pathways = (Pathway*) malloc(numReactions * sizeof (Pathway));
   for (int r = 0; r < numReactions; r++) {
      reactions.push_back(network.reactions[r]);
      pathways[r] = Pathway(r);
      reversible.push_back(network.reversible[r]);
      for (int m = 0; m < numMetabolites; m++) {
         pathways[r].setMetaboliteCoefficient(m, network.s[m][r]);
      }
   }
   for (int m = 0; m < numMetabolites; m++) {
      metabolites.push_back(network.metabolites[m]);
      externalMetabolites.push_back(network.external[m]);
   }
   for (unsigned int r = 0; r < network.reversiblePairs.size(); r++) {
      reversiblePairs.push_back(ReactionBitData());
      reversiblePairs[r].setBit(network.reversiblePairs[r], true);
      reversiblePairs[r].setBit(network.reversiblePairs[r] + 1, true);
   }
   EFMGenerator efmgenerator(pathways);
   efmgenerator.genenrateEFMs();
   free(pathways);
}

int main(int argc, char** argv) {
   if (argc == 2) {
      execute(argv[1]);
   } else {
      printf("Please specify the network file.");
   }
   return (EXIT_SUCCESS);
}

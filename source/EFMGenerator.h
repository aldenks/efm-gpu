#ifndef EFMGENERATOR_H
#define	EFMGENERATOR_H

#include <iostream>
#include "BitPatternTree.h"
#include "PathwayPool.h"
#include "ReversibleTree.h"
#include "Stopwatch.h"

using namespace std;

//Maximum number of combinations stored before processing of supersets.
#define MAX_COMBINATION_LIMIT 100000

class EFMGenerator {
private:
   //List of partial pathways
   PathwayPool pathways;
   Pathway* pathwaysPtr[MAX_COMBINATION_LIMIT];
   int combinationCount;
   //Total number of combinations generated
   long combinationsGenerated;
   //Total number combinations stored in the bitpattern tree
   //This is equal to total number of combinations - combinations rejected on the fly
   long combinationsStored;
   //Timer for removal of a metabolite
   Stopwatch timeMetaboliteRemvoed;
   BitPatternTree bptInput, bptOutput, bptNonpart;

   float format(double time);
   void removeUnusedInputs(int metaboliteIndex);
   void removeUnusedOutputs(int metaboliteIndex);
   void findNextMetaboliteToRemove(int& metaboliteIndex, int& inputCount, int& outputCount);
   void processExternals();
   void removeNextMetabolite();
   void printMetaboliteData();
   void reorderMetabolites(int index1, int index2);
   void generateCombinations(int metaboliteIndex, int inputCount, int outputCount);
   void generateCombinations(ReversibleTreeNode* input, ReversibleTreeNode* output);

public:
   EFMGenerator(Pathway* pths);
   ~EFMGenerator();
   void genenrateEFMs();
};

#endif	/* EFMGENERATOR_H */


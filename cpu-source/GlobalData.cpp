#include <cassert>
#include "BitPatternTreeNode.h"
#include "ReversibleTreeNode.h"

int numReactions;
int numMetabolites;
vector<bool> externalMetabolites;
vector<bool> reversible;
int numMetabolitesRemaining;
vector<string> reactions;
vector<string> metabolites;

int maxPoolSizeBPT;
BitPatternTreeNode* poolBPT;
int poolSizeBPT;

int maxPoolSizeRev;
ReversibleTreeNode* poolRev;
int poolSizeRev;

int reactionBitDataLength;
vector<ReactionBitData> reversiblePairs;
int BIT_COUNT_LOOKUP[] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

void initBPTNodePool() {
   maxPoolSizeBPT = MAX_BPT_TREE_NODES;
   poolBPT = (BitPatternTreeNode*) malloc(maxPoolSizeBPT * sizeof (BitPatternTreeNode));
   assert(poolBPT != 0);
   poolSizeBPT = 0;
}

void freeBPTNodePool() {
   free(poolBPT);
   poolBPT = 0;
}

void clearBPTNodePool() {
   poolSizeBPT = 0;
}

void* nextBPTNode() {
   poolSizeBPT++;
   assert(poolSizeBPT < maxPoolSizeBPT);
   poolBPT[poolSizeBPT - 1].init();
   return &poolBPT[poolSizeBPT - 1];
}

void initRevNodePool() {
   maxPoolSizeRev = MAX_REVERSIBLE_TREE_NODES;
   poolRev = (ReversibleTreeNode*) malloc(maxPoolSizeRev * sizeof (ReversibleTreeNode));
   assert(poolRev != 0);
   poolSizeRev = 0;
}

void freeRevNodePool() {
   free(poolRev);
   poolRev = 0;
}

void clearRevNodePool() {
   poolSizeRev = 0;
}

void* nextRevNode() {
   poolSizeRev++;
   assert(poolSizeRev < maxPoolSizeRev);
   return &poolRev[poolSizeRev - 1];
}

//Checks if a pathways contains reversible reaction pair

bool isValidPathway(ReactionBitData& rbd) {
   for (unsigned int i = 0; i < reversiblePairs.size(); i++) {
      if (rbd.isSuperSetOf(reversiblePairs[i])) {
         return false;
      }
   }
   return true;
}

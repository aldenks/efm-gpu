#ifndef REVERSIBLETREE_H
#define	REVERSIBLETREE_H

#include "ReversibleTreeNode.h"

using namespace std;

class ReversibleTree {
public:

   ReversibleTree(Pathway** list, int start, int end) {
      ReactionBitData bitsCannotbeUsed(true);
      getBitsCannotBeUsed(bitsCannotbeUsed, list, start, end);
      ReactionBitData bitsUsed;
      root = (ReversibleTreeNode*) nextRevNode();
      root->init(list, start, end, bitsCannotbeUsed, bitsUsed);
   }

   ~ReversibleTree() {
   }

   ReversibleTreeNode* getRoot() {
      return root;
   }


private:
   //Root of the tree
   ReversibleTreeNode* root;

   void getBitsCannotBeUsed(ReactionBitData& bitsCannotbeUsed, Pathway** list, int start, int end) {
      unsigned int count;
      for (int r = 0; r < numReactions; r++) {
         bitsCannotbeUsed.setBit(r, !reversible[r]);
         if (!bitsCannotbeUsed[r]) {
            count = 0;
            for (int p = start; p < end; p++) {
               if (list[p]->reactionBitData[r]) {
                  count++;
               } else {
                  count--;
               }
            }
            if (count < 0) {
               count = -count;
            }
         }
      }
   }
};

#endif	/* REVERSIBLETREE_H */


#ifndef BITPATTERNTREE_H
#define	BITPATTERNTREE_H

#include "BitPatternTreeNode.h"

class BitPatternTree {
public:
   BitPatternTreeNode* root;

   BitPatternTree() {
   }

   void init() {
      root = (BitPatternTreeNode*) nextBPTNode();
   }

   ~BitPatternTree() {
   }

   //Adds pathway to root node

   void addPathway(Pathway* p) {
      for (int i = 0; i < numReactions; i++) {
         bitsUsed[i] = false;
      }
      root->addPathway(p, bitsUsed);
   }

   //Checks if pathway is superset

   bool isSuperSet(ReactionBitData& rbd) {
      return root->isSuperSet(rbd);
   }

   bool isSuperSet(ReactionBitData& rbd, ReactionBitData& ref) {
      return root->isSuperSet(rbd, ref);
   }

   bool isSuperSet(Pathway* p) {
      return root->isSuperSet(p);
   }

   bool isSuperSet(Pathway* p, Pathway* ref) {
      return root->isSuperSet(p, ref);
   }
private:
   bool bitsUsed [MAX_REACTION_COUNT];
};

#endif	/* BITPATTERNTREE_H */


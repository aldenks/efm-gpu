#ifndef BITPATTERNTREENODE_H
#define	BITPATTERNTREENODE_H

#include "Pathway.h"

//Maximum number of pathways in leaf node
#define MAX_PATHWAYS  10

#define MAX_BPT_TREES 3
#define MAX_BPT_NODES_PER_TREE (MAX_REACTION_COUNT * MAX_REACTION_COUNT / MAX_PATHWAYS)
#define MAX_BPT_TREE_NODES (MAX_BPT_TREES * MAX_BPT_NODES_PER_TREE)

using namespace std;

void initBPTNodePool();
void freeBPTNodePool();
void* nextBPTNode();
void clearBPTNodePool();

class BitPatternTreeNode {
public:
   //List containing pathways
   Pathway* pathways[MAX_PATHWAYS + 1];

   BitPatternTreeNode() {
      init();
   }

   void init() {
      leaf = true;
      node0 = null;
      node1 = null;
      bit = -1;
      pathwayCount = 0;
      label.setAllBits();
   }

   ~BitPatternTreeNode() {
   }

   //Adds pathway

   void addPathway(Pathway* p, bool* bitsUsed) {
      //Check if node is leaf
      if (leaf) {
         pathways[pathwayCount++] = p;
         //Split if the number of pathways exceeds
         if (pathwayCount > MAX_PATHWAYS) {
            split(bitsUsed);
         }
         label.bitwiseAnd(p->reactionBitData);
      }//If node is not leaf
      else {
         bitsUsed[bit] = true;
         //If the bit is one store in node 1
         if (p->reactionBitData[bit]) {
            node1->addPathway(p, bitsUsed);
         } else {
            node0->addPathway(p, bitsUsed);
         }
         label.bitwiseAnd(p->reactionBitData);
      }
   }


   //Check if pathway is super set

   bool isSuperSet(ReactionBitData& rbd) {
      if (rbd.notAndEqualsNot(label)) {
         return false;
      }
      //Iterate through all the pathways in leaf node
      if (leaf) {
         bool subset, superset;
         for (unsigned int i = 0; i < pathwayCount; i++) {
            rbd.setOperations((pathways[i])->reactionBitData, subset, superset);
            //If set is p ignore it
            if (subset) {
               continue;
            }
            //Check if p is superset
            if (superset) {
               return true;
            }
         }
         return false;
      } else {
         //Select proper child node
         if (rbd[bit]) {
            if (node1->isSuperSet(rbd)) {
               return true;
            }
         }
         return node0->isSuperSet(rbd);
      }
   }
   //Check if pathway is super set

   bool isSuperSet(ReactionBitData& rbd, ReactionBitData& ref) {
      if (rbd.notAndEqualsNot(label)) {
         return false;
      }
      //Iterate through all the pathways in leaf node
      if (leaf) {
         for (unsigned int i = 0; i < pathwayCount; i++) {
            ReactionBitData pth = (pathways[i])->reactionBitData;
            //If set is ref ignore it
            if (pth == ref) {
               continue;
            }
            //Check if p is superset
            if (rbd.isSuperSetOf(pth)) {
               return true;
            }
         }
         return false;
      } else {
         //Select proper child node
         if (rbd[bit]) {
            if (node1->isSuperSet(rbd, ref)) {
               return true;
            }
         }
         return node0->isSuperSet(rbd, ref);
      }
   }

   //Check if pathway is super set

   bool isSuperSet(Pathway* p) {
      return isSuperSet(p->reactionBitData);
   }

   bool isSuperSet(Pathway* p, Pathway* ref) {
      return isSuperSet(p->reactionBitData, ref->reactionBitData);
   }

   BitPatternTreeNode* getNode0() {
      return node0;
   }

   BitPatternTreeNode* getNode1() {
      return node1;
   }

   bool isLeaf() {
      return leaf;
   }
private:
   //Current node is leaf if this is true
   bool leaf;
   //Node storing zero set
   BitPatternTreeNode* node0;
   //Node storing one set
   BitPatternTreeNode* node1;
   //Bit position used for decision
   int bit;
   //Label for the node
   ReactionBitData label;
   //Size of pathway vector
   unsigned int pathwayCount;

   //Split current node

   void split(bool* bitsUsed) {
      int splitCount [MAX_REACTION_COUNT];
      //Select the unused bit that has max number of ones and zeros
      for (int i = 0; i < numReactions; i++) {
         splitCount[i] = 0;
         if (bitsUsed[i]) {
            continue;
         }
         for (unsigned int j = 0; j < pathwayCount; j++) {
            if ((pathways[j])->reactionBitData[i]) {
               splitCount[i]++;
            } else {
               splitCount[i]--;
            }
         }
         if (splitCount[i] < 0) {
            splitCount[i] = -splitCount[i];
         }
      }
      int min = MAX_PATHWAYS;
      int height = 1;
      for (int i = 0; i < numReactions; i++) {
         if (bitsUsed[i]) {
            height++;
            continue;
         }
         if (min > splitCount[i]) {
            min = splitCount[i];
            bit = i;
         }
      }
      if (bit < 0) {
         return;
      }
      //Mark node as intermediate node
      leaf = false;
      node0 = (BitPatternTreeNode*) nextBPTNode();
      node1 = (BitPatternTreeNode*) nextBPTNode();
      //Copy pathways to new leaf nodes
      for (unsigned int i = 0; i < pathwayCount; i++) {
         if (pathways[i]->reactionBitData[bit]) {
            node1->pathways[node1->pathwayCount++] = pathways[i];
            node1->label.bitwiseAnd((pathways[i])->reactionBitData);
         } else {
            node0->pathways[node0->pathwayCount++] = pathways[i];
            node0->label.bitwiseAnd((pathways[i])->reactionBitData);
         }
      }
      pathwayCount = 0;
   }
};

#endif	/* BITPATTERNTREENODE_H */


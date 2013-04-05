#ifndef REVERSIBLETREENODE_H
#define	REVERSIBLETREENODE_H

using namespace std;

#define MAX_REVERSIBLE_TREES 2
#define MAX_REVERSIBLE_NODES_PER_TREE (MAX_REACTION_COUNT * 50)
#define MAX_REVERSIBLE_TREE_NODES (MAX_REVERSIBLE_TREES * MAX_REVERSIBLE_NODES_PER_TREE)

void initRevNodePool();
void freeRevNodePool();
void* nextRevNode();
void clearRevNodePool();

class ReversibleTreeNode {
public:

   ReversibleTreeNode() {
   }

   void init(Pathway** list, int s, int e, ReactionBitData& bitsCannotBeUsed, ReactionBitData& bUsed) {
      pathways = list;
      leaf = true;
      start = s;
      end = e;
      node0 = node1 = null;
      bitsUsed = bUsed;
      //Check if the all the pathways in this partition are valid
      if (!isValidPathway(bitsUsed)) {
         end = start;
         return;
      }
      bit = getSplitBit(bitsCannotBeUsed);
      //Check if further partitioning is required
      if (bit < 0) {
         checkAllUnusedBits(bitsCannotBeUsed);
         return;
      }
      //Split the data
      ReactionBitData updatedBitsCannotBeUsed(bitsCannotBeUsed);
      updatedBitsCannotBeUsed.setBit(bit, true);
      leaf = false;
      int middle = partition();
      ReactionBitData bUsed0(bUsed);
      node0 = (ReversibleTreeNode*) nextRevNode();
      node0->init(pathways, start, middle, updatedBitsCannotBeUsed, bUsed0);
      bUsed.setBit(bit, true);
      node1 = (ReversibleTreeNode*) nextRevNode();
      node1->init(pathways, middle, end, updatedBitsCannotBeUsed, bUsed);
   }

   ~ReversibleTreeNode() {
   }

   int size() {
      return end - start;
   }

   bool isLeaf() {
      return leaf;
   }

   ReversibleTreeNode* getNode0() {
      return node0;
   }

   ReversibleTreeNode* getNode1() {
      return node1;
   }

   ReactionBitData& getBitsUsed() {
      return bitsUsed;
   }

   int getStart() {
      return start;
   }

   int getEnd() {
      return end;
   }

   Pathway* getPathway(int p) {
      return pathways[p];
   }



private:
   //Pathway list
   Pathway** pathways;
   //Flag indicating if the current node is a leaf
   bool leaf;
   //Starting and ending index of pathways in the list
   int start, end;
   //Children nodes if current pathway is not a leaf
   //node0 represents the subtree in which the reaction used as split point is inactive
   //node1 represents the subtree in which the reaction used as split point is active
   ReversibleTreeNode* node0;
   ReversibleTreeNode* node1;
   //Reaction used as split point
   int bit;
   //Reactions that have been used for splitting the tree so far from root of the tree to this node
   ReactionBitData bitsUsed;

   void checkAllUnusedBits(ReactionBitData& bitsCannotBeUsed) {
      for (int r = 0; r < numReactions; r++) {
         if (!bitsCannotBeUsed[r]) {
            if (pathways[start]->reactionBitData[r]) {
               bitsUsed.setBit(r, true);
            }
         }
      }
   }
   //Finds reaction that will be used to split pathways at current node

   int getSplitBit(ReactionBitData& bitsCannotBeUsed) {
      int bit = -1;
      int count;
      int minCount = size();
      for (int r = 0; r < numReactions; r++) {
         if (bitsCannotBeUsed[r]) {
            continue;
         }
         count = 0;
         for (int p = start; p < end; p++) {
            if (pathways[p]->reactionBitData[r]) {
               count++;
            } else {
               count--;
            }
         }
         if (count < 0) {
            count = -count;
         }
         if (count < minCount) {
            minCount = count;
            bit = r;
         }
      }
      return bit;
   }
   //Sorts pathways such that pathways in which the split reaction is inactive appears before the pathways in which the split reaction is active

   int partition() {
      int p0, p1;
      Pathway* p;
      for (p0 = start, p1 = end - 1; p0 <= p1; p0++, p1--) {
         //Pathway p0 has bit set
         if (pathways[p0]->reactionBitData[bit]) {
            //Pathway p1 has bit set
            if (pathways[p1]->reactionBitData[bit]) {
               //Do not skip p0
               p0--;
            }//Pathway p1 has bit reset
            else {
               p = pathways[p0];
               pathways[p0] = pathways[p1];
               pathways[p1] = p;
            }
         }//Pathway p0 has bit reset
         else {
            //Pathway p1 has bit reset
            if (!(pathways[p1]->reactionBitData[bit])) {
               //Do not skip p1
               p1++;
            }
         }
      }
      return p0;
   }
};

#endif	/* REVERSIBLETREENODE_H */


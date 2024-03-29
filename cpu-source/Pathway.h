#ifndef PATHWAY_H
#define	PATHWAY_H

#define MAX_METABOLITES 100

#include "ReactionBitData.h"

class Pathway {
private:
   //Parent pathways that were combined to generate this pathway
   Pathway* parent1;
   Pathway* parent2;
public:
   //Represents active reactions in the pathway
   ReactionBitData reactionBitData;
   //Metabolite coefficients for the pathways
   double metaboliteCoefficients [MAX_METABOLITES];
   int metaboliteCount;

   Pathway() {
      parent1 = null;
      parent2 = null;
      metaboliteCount = 0;
   }

   Pathway(int reaction) {
      parent1 = null;
      parent2 = null;
      metaboliteCount = numMetabolites;
      for (int i = 0; i < metaboliteCount; i++) {
         metaboliteCoefficients[i] = 0;
      }
      reactionBitData.setBit(reaction, true);
   }

   ~Pathway() {
   }

   void setParents(Pathway* p1, Pathway* p2) {
      parent1 = p1;
      parent2 = p2;
      metaboliteCount = 0;
      rbdOr(reactionBitData, p1->reactionBitData, p2->reactionBitData);
   }

   void setMetaboliteCoefficient(int m, double c) {
      metaboliteCoefficients[m] = c;
   }

   void updateMetaboliteCoefficients(int m) {
      metaboliteCount = numMetabolitesRemaining;
      double scale = -parent1->metaboliteCoefficients[m] / parent2->metaboliteCoefficients[m];
      //This is done to keep co-efficients greater than 1
      if (scale >= 1) {
         for (int i = 0; i < metaboliteCount; i++) {
            metaboliteCoefficients[i] = parent1->metaboliteCoefficients[i] + scale * parent2->metaboliteCoefficients[i];
         }
      } else {
         scale = 1 / scale;
         for (int i = 0; i < metaboliteCount; i++) {
            metaboliteCoefficients[i] = scale * parent1->metaboliteCoefficients[i] + parent2->metaboliteCoefficients[i];
         }
      }
      parent1 = parent2 = null;
   }

   bool isInput(int m) {
      return (metaboliteCoefficients[m] < NEG_ZERO);
   }

   bool isOutput(int m) {
      return (metaboliteCoefficients[m] > ZERO);
   }

   bool isSupersetOf(Pathway& p) {
      return reactionBitData.isSuperSetOf(p.reactionBitData);
   }
};

bool isValidPathway(ReactionBitData& rbd);

#endif	/* PATHWAY_H */

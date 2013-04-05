#ifndef REACTIONBITDATA_H
#define	REACTIONBITDATA_H A

#include "GlobalData.h"

#define MAX_REACTION_COUNT 320
//Length of reaction bit data array
extern int reactionBitDataLength;
extern int BIT_COUNT_LOOKUP[];

#define CARDINALITY(n) BIT_COUNT_LOOKUP[n & 0xFF] + BIT_COUNT_LOOKUP[(n >> 8) & 0xFF] + BIT_COUNT_LOOKUP[(n >> 16) & 0xFF] + BIT_COUNT_LOOKUP[(n >> 24) & 0xFF]

class ReactionBitData {
public:
   unsigned int data [MAX_REACTION_COUNT / sizeof (int)];

   ReactionBitData() {
      for (int i = 0; i < reactionBitDataLength; i++) {
         data[i] = 0;
      }
   }

   ReactionBitData(bool allbits) {
      unsigned int value = allbits ? -1 : 0;
      for (int i = 0; i < reactionBitDataLength; i++) {
         data[i] = value;
      }
   }

   ReactionBitData(const ReactionBitData& bd) {
      for (int i = 0; i < reactionBitDataLength; i++) {
         data[i] = bd.data[i];
      }
   }

   ~ReactionBitData() {
   }

   void setBit(int bit, bool value) {
      int index = bit >> 5;
      bit &= 0x1F;
      unsigned int mask = 1 << bit;
      if (value) {
         data[index] |= mask;
      } else {
         data[index] &= (mask ^ -1);
      }
   }

   void setAllBits() {
      for (int i = 0; i < reactionBitDataLength; i++) {
         data[i] = -1;
      }
   }

   const bool operator[] (const int bit) const {
      return (data[bit >> 5] & (1 << (bit & 0x1F)));
   }

   void bitwiseAnd(const ReactionBitData& obj) {
      for (int i = 0; i < reactionBitDataLength; i++) {
         data[i] = data[i] & obj.data[i];
      }
   }

   bool isSuperSetOf(const ReactionBitData& obj) {
      for (int i = 0; i < reactionBitDataLength; i++) {
         if ((data[i] & obj.data[i]) != obj.data[i]) {
            return false;
         }
      }
      return true;
   }

   void setOperations(const ReactionBitData& obj, bool &subset, bool &superset) {
      unsigned int andVal;
      subset = true;
      superset = true;
      for (int i = 0; (i < reactionBitDataLength) && (subset || superset); i++) {
         andVal = data[i] & obj.data[i];
         if (andVal != data[i]) {
            subset = false;
         }
         if (andVal != obj.data[i]) {
            superset = false;
         }
      }
   }

   bool notAndEqualsNot(const ReactionBitData& rbd) {
      for (int i = 0; i < reactionBitDataLength; i++) {
         if (((-1 - data[i]) & rbd.data[i]) > 0) {
            return true;
         }
      }
      return false;
   }

   bool operator==(const ReactionBitData& obj) {
      for (int i = 0; i < reactionBitDataLength; i++) {
         if (data[i] != obj.data[i]) {
            return false;
         }
      }
      return true;
   }

   ReactionBitData& operator=(const ReactionBitData& bd) {
      for (int i = 0; i < reactionBitDataLength; i++) {
         data[i] = bd.data[i];
      }
      return *this;
   }

   int getCardinality() {
      int c = 0;
      for (int i = 0; i < reactionBitDataLength; i++) {
         c += CARDINALITY(data[i]);
      }
      return c;
   }

private:
};

inline void rbdAnd(ReactionBitData& output, const ReactionBitData& a, const ReactionBitData& b) {
   for (int i = 0; i < reactionBitDataLength; i++) {
      output.data[i] = a.data[i] & b.data[i];
   }
}

inline void rbdOr(ReactionBitData& output, const ReactionBitData& a, const ReactionBitData& b) {
   for (int i = 0; i < reactionBitDataLength; i++) {
      output.data[i] = a.data[i] | b.data[i];
   }
}

//Reversible reaction pairs
extern vector<ReactionBitData> reversiblePairs;

#endif	/* REACTIONBITDATA_H */

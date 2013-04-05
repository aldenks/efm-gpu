#ifndef PATHWAYPOOL_H
#define	PATHWAYPOOL_H

#include <cassert>

#define MAX_PATHWAY_POOL_SIZE 1000000

class PathwayPool {
public:

   PathwayPool() {
      pathways = (Pathway*) malloc(MAX_PATHWAY_POOL_SIZE * sizeof (Pathway));
      assert(pathways != 0);
      nCount = 0;
   }

   ~PathwayPool() {
      free(pathways);
      pathways = 0;
   }

   unsigned int size() {
      return nCount;
   }

   Pathway& operator[] (const unsigned int nIndex) {
      assert(nIndex >= 0 && nIndex < nCount);
      return pathways[nIndex];
   }

   void operator<<(Pathway& pathway) {
      pathways[nCount] = pathway;
      nCount++;
   }

   Pathway* operator ++() {
      return &pathways[nCount++];
   }

   void operator --() {
      nCount--;
   }

   void remove(const unsigned int nIndex) {
      nCount--;
      pathways[nIndex] = pathways[nCount];
   }

private:
   Pathway* pathways;
   unsigned int nCount;
};

#endif	/* PATHWAYPOOL_H */

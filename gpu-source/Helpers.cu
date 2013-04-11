#include "Helpers.h"
#include "assert.h"

void markMetaboliteBalanced(int metabolite, bool* d_balancedMetabolites) {
   // hack which relies on this assertion being true. should compile away to nothing.
   assert(sizeof(bool) == 1); 
   cudaMemset(d_balancedMetabolites + metabolite, true, 1); 
}

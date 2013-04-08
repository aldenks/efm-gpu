#include "EFMGenerator.h"

//Generates EFMs

void generateEFMs() {
   while (remainingMetabolites > 0) {
      //Call a kernel to report the index of metabolite to be removed
      //Copy the bit vectors from gpu to cpu
      //Sort the bit vectors for inputs outputs and non-participating pathways
      //Copy the bit vectors from cpu to gpu
      //Setup bins
      //Call kernel to generate combinations
      remainingMetabolites--;
   }
}
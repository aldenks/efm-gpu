#include "GenerateCombos.h"

__device__
void scalePathway(float* result, float* p1, float* p2, float scale, int metaboliteCount){
   for(int i = 0; i < metaboliteCount; i++){
      result[i] = p1[i] + p2[i] * scale;
   }
}

__global__
void generateCombinations(int* bins, int* indicies, int inputIndex, int numberOfBins, int metabolite, int metaboliteCount, BinaryVector *reactions, float *metaboliteCoefficients){
   int tid = blockIdx.x + blockDim.x + threadIdx.x;
   if(tid >= numberOfBins){
      return;
   }
   int writeIndex = indicies[tid];
   bins += tid;
   int writeCount = *bins;
   inputIndex = circularIndex(inputIndex + tid);
   BinaryVector inputReaction = reactions[inputIndex];
   int outputIndex;
   float met1, met2;
   for(int i = 0; i < writeCount; i++, writeIndex = circularIndex(writeIndex + 1)){
      bins += numberOfBins;
      outputIndex = *bins;
      reactions[writeIndex] = inputReaction | reactions[outputIndex];
      met1 = metaboliteCoefficients[metaboliteCount * inputIndex + metabolite];
      met2 = metaboliteCoefficients[metaboliteCount * outputIndex + metabolite];
      if(met1 < met2){
         scalePathway(metaboliteCoefficients + metaboliteCount * writeIndex, metaboliteCoefficients + metaboliteCount * inputIndex, metaboliteCoefficients + metaboliteCount * outputIndex, met1/met2, metaboliteCount);
      }else{
         scalePathway(metaboliteCoefficients + metaboliteCount * writeIndex, metaboliteCoefficients + metaboliteCount * outputIndex, metaboliteCoefficients + metaboliteCount * inputIndex, met2/met1, metaboliteCount);
      }
   }
}

int generateCombinations(int metabolite, int numberOfBins, int nextFreePathwayIndex){
   cudaMemcpy(h_combinationBinCounts, d_combinationBins, numberOfBins * sizeof(int), cudaMemcpyDeviceToHost);
   h_combinationWriteIndices[0] = nextFreePathwayIndex;
   int newComboCount = h_combinationBinCounts[0];
   printf("%i -- %i combos\n", 0, h_combinationBinCounts[0]);
   for(int i = 1; i < numberOfBins; i++){
      printf("%i -- %i combos\n", i, h_combinationBinCounts[i]);
      h_combinationWriteIndices[i] = h_combinationWriteIndices[i - 1] + h_combinationBinCounts[i - 1];
      newComboCount += h_combinationBinCounts[i];
   }
   cudaMemcpy(d_combinationWriteIndices, h_combinationWriteIndices, numberOfBins * sizeof(int), cudaMemcpyHostToDevice);
   int numBlocks = (numberOfBins / MAX_THREADS_PER_BLOCK) + 1;
   generateCombinations << < numBlocks, MAX_THREADS_PER_BLOCK >> > (d_combinationBins, d_combinationWriteIndices, pathwayStartIndex, numberOfBins, metabolite, metaboliteCount, d_binaryVectors, d_metaboliteCoefficients);
   return newComboCount;
}

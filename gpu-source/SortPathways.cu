#include "SortPathways.h"

__global__
void sortInputPathways(BinaryVector *reactions, float *metaboliteCoefficients, int startIndex, int numToSort, int numReactions, int numberOfMetabolites, int metaboliteToRemove) {

   int tid = blockIdx.x + blockDim.x + threadIdx.x;

   if (tid == 0) {
      int start = startIndex;
      int end = startIndex + numToSort;

      //While the pointers do not overlap
      while (start < end) {
	bool is_input_1 = metaboliteCoefficients[circularIndex(start) * numberOfMetabolites + metaboliteToRemove] < NEG_ZERO;
	bool is_input_2 = metaboliteCoefficients[circularIndex(end) * numberOfMetabolites + metaboliteToRemove] < NEG_ZERO;
         if (is_input_1) {
            //Skip this one
            start++;
         } else if (is_input_2) {
            //swap the two reactions
	   BinaryVector temp = reactions[circularIndex(end)];
	   reactions[circularIndex(end)] = reactions[circularIndex(start)];
	   reactions[circularIndex(start)] = temp;

            //swap the two metaboliteCoefficients for pathways
            float tempCoefficient;
            for (int i = 0; i < numberOfMetabolites; ++i) {
	      tempCoefficient = metaboliteCoefficients[circularIndex(end) * numberOfMetabolites + i];
	      metaboliteCoefficients[circularIndex(end) * numberOfMetabolites + i] = metaboliteCoefficients[circularIndex(start) * numberOfMetabolites + i];
	      metaboliteCoefficients[circularIndex(start) * numberOfMetabolites + i] = tempCoefficient;
            }
            //move forward
            start++;
            end--;
         } else {
            //Not an input, don't care
            end--;
         }
      }
   }
}

__global__
void sortOutputPathways(BinaryVector *reactions, float *metaboliteCoefficients, int startIndex, int numToSort, int numReactions, int numberOfMetabolites, int metaboliteToRemove) {
   int tid = blockIdx.x + blockDim.x + threadIdx.x;
   if (tid == 0) {
      int start = startIndex;
      int end = startIndex + numToSort;
      int startIndex = circularIndex(start);
      int endIndex = circularIndex(end);
      //While the pointers do not overlap
      while (start < end) {
         bool is_output_1 = metaboliteCoefficients[startIndex * numberOfMetabolites +
                 metaboliteToRemove] > ZERO;
         bool is_output_2 = metaboliteCoefficients[endIndex * numberOfMetabolites +
                 metaboliteToRemove] > ZERO;
         if (is_output_1) {
            //Skip this one
            start++;
                        startIndex = circularIndex(start);
         }
         if (is_output_2) {
            //swap the two reactions
            BinaryVector temp = reactions[endIndex];
            reactions[endIndex] = reactions[startIndex];
            reactions[startIndex] = temp;

            //swap the two metaboliteCoefficients for pathways
            float tempCoefficient;
            for (int i = 0; i < numberOfMetabolites; ++i) {
               tempCoefficient = metaboliteCoefficients[endIndex * numberOfMetabolites + i];
               metaboliteCoefficients[endIndex * numberOfMetabolites + i] =
                       metaboliteCoefficients[startIndex * numberOfMetabolites + i];
               metaboliteCoefficients[startIndex * numberOfMetabolites + i] = tempCoefficient;
            }
            //move forward
            start++;
            end--;

            startIndex = circularIndex(start);
            endIndex = circularIndex(end);
         } else {
            //Not an output, don't care
            end--;
                        endIndex = circularIndex(end);
         }
      }
   }
}

__global__
void dependencyCheck(BinaryVector *reactions, int *bins, int batch_size, int num_inputs, int output_start, int non_part_start, int pathwayCounts, int pathwayStartIndex) {
   int tid = blockIdx.x * blockDim.x + threadIdx.x;
   int count = 0;
   if (tid < num_inputs) {
      bins[tid] = 0;
      BinaryVector input = reactions[circularIndex(pathwayStartIndex + tid)];
      BinaryVector output, combo, pathway;
      bool is_unique_and_independent = true;
      for (int i = 0; i < batch_size && output_start + i < non_part_start; ++i) {
         output = reactions[circularIndex(output_start + i)];
         combo = input | output;
         for (int j = 0; is_unique_and_independent && j < pathwayCounts; ++j) {
            if (j == tid) {
               continue; //skip this input
            }
            if (j == output_start + i) {
               continue; //skip this output
            }
            pathway = reactions[circularIndex(j)];
            if (pathway == combo) {
               //TODO: how can we prevent duplicates?
               is_unique_and_independent = false;
               break;
            }
            is_unique_and_independent = ((combo & pathway) != pathway);
         }
         if (is_unique_and_independent) {
            count++;
            bins[num_inputs * (count) + tid] = i;
            
         }
      }
      /*
      //Assumes first item is 0 to start
      if (circularIndex(output_start) == circularIndex(pathwayStartIndex + num_inputs)) {
         bins[tid] = 0;
      }*/
      bins[tid] = count;
   }
}

__global__
void checkSort(BinaryVector *reactions, float *metaboliteCoefficients, int bufferStart, int numPathways, int numMetabolites, int metaboliteToRemove, int *sorts) {
   int tid = blockIdx.x * blockDim.x + threadIdx.x;
   if (tid < numPathways) {
      int startIndex = circularIndex(bufferStart + tid);
      if (metaboliteCoefficients[startIndex * numMetabolites + metaboliteToRemove] < NEG_ZERO) {
         sorts[startIndex] = 0;
      } else if (metaboliteCoefficients[startIndex * numMetabolites + metaboliteToRemove] > ZERO) {
         sorts[startIndex] = 1;
      } else {
         sorts[startIndex] = 2;
      }
   }
}

void sortInputsOutputs(int pathwayCounts, int metaboliteCount, int numInputs, int numOutputs, int metaboliteToRemove) {
   //call the kernel on inputs and outputs
   printf("sortInputsOutputs: pathwayCount=%d, metaboliteCount=%d, numInputs=%d, numbOutputs=%d, metabolite=%d\n",pathwayCounts, metaboliteCount, numInputs, numOutputs, metaboliteToRemove);
   int numBlocks = (pathwayCounts / MAX_THREADS_PER_BLOCK) + 1;
   sortInputPathways << < numBlocks, MAX_THREADS_PER_BLOCK >> > (d_binaryVectors, d_metaboliteCoefficients, pathwayStartIndex, pathwayCounts, pathwayCounts, metaboliteCount, metaboliteToRemove);

   numBlocks = ((pathwayCounts - numInputs) / MAX_THREADS_PER_BLOCK) + 1;
   sortOutputPathways << < numBlocks, MAX_THREADS_PER_BLOCK >> > (d_binaryVectors, d_metaboliteCoefficients, pathwayStartIndex + numInputs, pathwayCounts - numInputs, pathwayCounts, metaboliteCount, metaboliteToRemove);

   int *h_sorts = (int *) malloc(sizeof (int) * MAX_PATHWAYS);
   for (int i = 0; i < MAX_PATHWAYS; ++i) {
      h_sorts[i] = -1;
   }

   int *d_sorts = NULL;
   cudaMalloc((void **) &d_sorts, MAX_PATHWAYS * sizeof (int));

   cudaError error;
   error = cudaMemcpy(d_sorts, h_sorts, sizeof (int)*MAX_PATHWAYS, cudaMemcpyHostToDevice);
   if (error != cudaSuccess)
      fprintf(stderr, "Error in copying for check sort HTD\n");
   numBlocks = (pathwayCount / MAX_THREADS_PER_BLOCK) + 1;
   checkSort << < numBlocks, MAX_THREADS_PER_BLOCK >> > (d_binaryVectors, d_metaboliteCoefficients, pathwayStartIndex, pathwayCount, metaboliteCount, metaboliteToRemove, d_sorts);

   error = cudaMemcpy(h_sorts, d_sorts, sizeof (int)*MAX_PATHWAYS, cudaMemcpyDeviceToHost);
   if (error != cudaSuccess)
      fprintf(stderr, "Error in copying for check sort DTH\n");

   for (int i = pathwayStartIndex; i < pathwayStartIndex + pathwayCount; i++) {
      fprintf(stderr, "Pathway %i is %i\n", circularIndex(i), h_sorts[circularIndex(i)]);
   }
/*
   for (int i = 0; i < MAX_PATHWAYS; ++i) {
      if (h_sorts[i] > -1) {
         fprintf(stderr, "Found something.");
      }
   }
*/
   cudaFree(d_sorts);
   free(h_sorts);
}

void dependencyCheck(int numInputs, int numOutputs, int batch_number) {
   int numBlocks = (numInputs / MAX_THREADS_PER_BLOCK) + 1;
   printf("dependencyCheck: numInputs=%d, numOutputs=%d, batch_number=%d, numBlocks=%d, batchSize=%d\n", numInputs, numOutputs, batch_number, numBlocks, batchSize);
   dependencyCheck << < numBlocks, MAX_THREADS_PER_BLOCK >> > (d_binaryVectors, d_combinationBins, batchSize, numInputs,
           pathwayStartIndex + numInputs + (batch_number * batchSize), //start of next batch of outputs
           pathwayStartIndex + numInputs + numOutputs, //start of non-participating
           pathwayCount, pathwayStartIndex);
}

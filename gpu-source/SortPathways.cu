#include "SortPathways.h"

__global__
void sortInputPathways(BinaryVector *reactions, float *metaboliteCoefficients, int input_start, int numToSort, int numberOfMetabolites, int metaboliteToRemove) {

   int tid = blockIdx.x * blockDim.x + threadIdx.x;

   if (tid == 0) {
      int start = input_start;

      //Need to remove 1 so that end represents the last index
      int end = input_start + numToSort - 1;

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
void sortOutputPathways(BinaryVector *reactions, float *metaboliteCoefficients, int output_start, int numToSort, int numberOfMetabolites, int metaboliteToRemove) {
   int tid = blockIdx.x * blockDim.x + threadIdx.x;

   if (tid == 0) {
      int start = output_start;

      //Need to remove 1 so that end represents the last index
      int end = output_start + numToSort - 1;

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
	      is_unique_and_independent = false; //Found a duplicate
               break;
            }

            is_unique_and_independent = ((combo & pathway) != pathway);
         }

         if (is_unique_and_independent) {
            count++;
            bins[num_inputs * (count) + tid] = circularIndex(output_start + i);
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
void checkSort(float *metaboliteCoefficients, int pathwayStartIndex, int numPathways, int numMetabolites, int metaboliteToRemove, int *sorts) {
   int tid = blockIdx.x * blockDim.x + threadIdx.x;
   //if (tid < numPathways) {
   if(tid == 0){
     //int startIndex = circularIndex(pathwayStartIndex + tid);
     for(int i = 0; i < numPathways; ++i){
       int startIndex = circularIndex(pathwayStartIndex + i);
       if (metaboliteCoefficients[startIndex * numMetabolites + metaboliteToRemove] < NEG_ZERO) {
         sorts[startIndex] = -1;
       } else if (metaboliteCoefficients[startIndex * numMetabolites + metaboliteToRemove] > ZERO) {
         sorts[startIndex] = 1;
       } else {
         sorts[startIndex] = 0;
       }
     }
   }
}

void checkSorting(int *h_sorts, int *d_sorts, int metaboliteToRemove){
   for (int i = 0; i < MAX_PATHWAYS; ++i) {
      h_sorts[i] = -2;
   }

   cudaError error;
   error = cudaMemcpy(d_sorts, h_sorts, sizeof (int)*MAX_PATHWAYS, cudaMemcpyHostToDevice);
   if (error != cudaSuccess)
      fprintf(stderr, "Error in copying for check sort HTD\n");

   int numBlocks = (pathwayCount / MAX_THREADS_PER_BLOCK) + 1;
   checkSort << < numBlocks, MAX_THREADS_PER_BLOCK >> > (d_metaboliteCoefficients, pathwayStartIndex, pathwayCount, metaboliteCount, metaboliteToRemove, d_sorts);

   error = cudaMemcpy(h_sorts, d_sorts, sizeof (int)*MAX_PATHWAYS, cudaMemcpyDeviceToHost);
   if (error != cudaSuccess)
      fprintf(stderr, "Error in copying for check sort DTH\n");

   for (int i = pathwayStartIndex; i < pathwayStartIndex + pathwayCount; i++) {
      fprintf(stderr, "Pathway %i is %i\n", circularIndex(i), h_sorts[circularIndex(i)]);
   }
}

void sortInputsOutputs(int numInputs, int numOutputs, int metaboliteToRemove) {
   //call the kernel on inputs and outputs
   printf("sortInputsOutputs: pathwayCount=%d, metaboliteCount=%d, numInputs=%d, numOutputs=%d, metabolite=%d\n",pathwayCount, metaboliteCount, numInputs, numOutputs, metaboliteToRemove);
   int numBlocks = (pathwayCount / MAX_THREADS_PER_BLOCK) + 1;

   int *h_sorts = (int *) malloc(sizeof (int) * MAX_PATHWAYS);
   int *d_sorts = NULL;
   cudaMalloc((void **) &d_sorts, MAX_PATHWAYS * sizeof (int));

   fprintf(stderr, "Before sort\n");
   checkSorting(h_sorts, d_sorts, metaboliteToRemove);

   sortInputPathways << < 1, 32 >> > (d_binaryVectors, d_metaboliteCoefficients, pathwayStartIndex, pathwayCount, metaboliteCount, metaboliteToRemove);

   fprintf(stderr, "Inputs sorted\n");
   checkSorting(h_sorts, d_sorts, metaboliteToRemove);

   numBlocks = ((pathwayCount - numInputs) / MAX_THREADS_PER_BLOCK) + 1;
   sortOutputPathways << < 1, 32 >> > (d_binaryVectors, d_metaboliteCoefficients, pathwayStartIndex + numInputs, pathwayCount - numInputs, metaboliteCount, metaboliteToRemove);

   fprintf(stderr, "Inputs and outputs sorted\n");
   checkSorting(h_sorts, d_sorts, metaboliteToRemove);

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

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

	//Recall that a metabolite is an input to a reaction if its input is negative
	bool is_input_1 = metaboliteCoefficients[circularIndex(start) * numberOfMetabolites + metaboliteToRemove] < NEG_ZERO;
	bool is_input_2 = metaboliteCoefficients[circularIndex(end) * numberOfMetabolites + metaboliteToRemove] < NEG_ZERO;

         if (is_input_1) {
            //Skip this one since it is in the right bucket
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
            //Not an input, already in the right bucket
            end--;
         }
      }
   }
}

//Assumes that sortInputPathways has already been called
__global__
void sortOutputPathways(BinaryVector *reactions, float *metaboliteCoefficients, 
			int output_start, int numToSort, int numberOfMetabolites, 
			int metaboliteToRemove) {
   int tid = blockIdx.x * blockDim.x + threadIdx.x;

   if (tid == 0) {
      int start = output_start;

      //Need to remove 1 so that end represents the last index
      int end = output_start + numToSort - 1;

      int startIndex = circularIndex(start);
      int endIndex = circularIndex(end);

      //While the pointers do not overlap
      while (start < end) {
	//Recall that a metabolite is an input to a reaction if its input is positive
         bool is_output_1 = metaboliteCoefficients[startIndex * numberOfMetabolites +
                 metaboliteToRemove] > ZERO;
         bool is_output_2 = metaboliteCoefficients[endIndex * numberOfMetabolites +
                 metaboliteToRemove] > ZERO;
         if (is_output_1) {
            //Skip this one since it is already sorted
            start++;
            startIndex = circularIndex(start);
         } else if (is_output_2) {
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
            //Not an output, so already sorted
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

      //Need to calculate start of this batch of outputs in circular buffer 
      int buffer_output_start = circularIndex(pathwayStartIndex + output_start);

      //Keep generating combos until either
      //  1) batch is finished OR
      //  2) reached the non participating pathways (end of output)
      for (int i = 0; i < batch_size && (output_start + i) < non_part_start; ++i) {
         output = reactions[circularIndex(buffer_output_start + i)];
         combo = input | output;

	 // We are trying to prove that combo is NOT independent
         bool is_unique_and_independent = true;

         for (int j = 0; is_unique_and_independent && j < pathwayCounts; ++j) {
	   //If INDEX of this pathway is the same as INDEX of input, skip it
            if (j == tid) {
               continue; 
            }

	   //If INDEX of this pathway is the same as INDEX of output, skip it
            if (j == output_start + i) {
               continue;
            }

            pathway = reactions[circularIndex(pathwayStartIndex + j)];

            if (pathway == combo) {
	      //Found a duplicate.
	      //NOTE: The pathway that is already in the array is considered
	      //      the "original" in this case.
               is_unique_and_independent = false; 
               break;
            }
	    
	    //Check if combo is not dependent on pathway
            is_unique_and_independent = ((combo & pathway) != pathway);
         }

	 //If unique and independent, store the output index into the bin
         if (is_unique_and_independent) {
            count++;
            bins[num_inputs * (count) + tid] = circularIndex(buffer_output_start + i);
         }
      }

      bins[tid] = count;
   }
}

// Takes the metabolite coefficients, number of pathways to consider, starting index in
// the circular buffer, number of metabolites, the metabolite we are considering, and
// an array of ints.

// Writes the following to the sorts array in the index corresponding to a pathway in
// metaboliteCoefficients:
//    - -1 if the metabolite is an input to that pathway
//    - 1 if the metabolite is an output to that pathway
//    - 0 otherwise (non-participating)
__global__
void checkSort(float *metaboliteCoefficients, int pathwayStartIndex, int numPathways, 
	       int numMetabolites, int metaboliteToRemove, int *sorts) {
   int tid = blockIdx.x * blockDim.x + threadIdx.x;
   if(tid == 0){
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

//Takes the metabolite we wish to check

//Prints every pathway considered in this iteration and whether they are
//an input, output, or non-participating.

//NOTE: This is used to ensure the buckets were properly made after sorting.
void checkSorting(int metaboliteToRemove){
   int *h_sorts = (int *) malloc(sizeof (int) * MAX_PATHWAYS);

   int *d_sorts = NULL;
   cudaError error;
   error = cudaMalloc((void **) &d_sorts, MAX_PATHWAYS * sizeof (int));
   if (error != cudaSuccess) {
      fprintf(stderr, "SortPathways.cu:checkSorting() Unable to allocate memory for d_sorts\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));

      free(h_sorts);
      return;
   }

   //Given default input to allow others to do some other checks
   for (int i = 0; i < MAX_PATHWAYS; ++i) {
      h_sorts[i] = -2;
   }


   error = cudaMemcpy(d_sorts, h_sorts, sizeof (int)*MAX_PATHWAYS, cudaMemcpyHostToDevice);
   if (error != cudaSuccess) {
      fprintf(stderr, "SortPathways.cu:checkSorting() Unable to copy h_sorts over to GPU\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));

      cudaFree(d_sorts);
      free(h_sorts);
      return;
   }

   // Dispatched for greater parallel exploitation, but kernel is currently a single thread
   int numBlocks = (pathwayCount / MAX_THREADS_PER_BLOCK) + 1;
   checkSort << < numBlocks, MAX_THREADS_PER_BLOCK >> > (d_metaboliteCoefficients, pathwayStartIndex, pathwayCount, metaboliteCount, metaboliteToRemove, d_sorts);

   error = cudaMemcpy(h_sorts, d_sorts, sizeof (int)*MAX_PATHWAYS, cudaMemcpyDeviceToHost);
   if (error != cudaSuccess) {
      fprintf(stderr, "SortPathways.cu:checkSorting() Unable to copy d_sorts over to CPU\n");
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(error));

      cudaFree(d_sorts);
      free(h_sorts);
      return;
   }

   for (int i = pathwayStartIndex; i < pathwayStartIndex + pathwayCount; i++) {
      fprintf(stderr, "Pathway %i is %i\n", circularIndex(i), h_sorts[circularIndex(i)]);
   }

   cudaFree(d_sorts);
   free(h_sorts);
}

// Calls the two sorting kernels, sortInputPathways() and sortOutputPathways()
void sortInputsOutputs(int numInputs, int numOutputs, int metaboliteToRemove) {
   //printf("sortInputsOutputs: pathwayCount=%d, metaboliteCount=%d, numInputs=%d, numOutputs=%d, metabolite=%d\n",pathwayCount, metaboliteCount, numInputs, numOutputs, metaboliteToRemove);

  // Sort all pathways into input bucket and other bucket.
  // NOTE: We only use one thread, so we dedicate one warp to it.
   sortInputPathways << < 1, 32 >> > (d_binaryVectors, d_metaboliteCoefficients, 
				      pathwayStartIndex, //Index to start sorting
				      pathwayCount,      //Number of pathways to sort
				      metaboliteCount, metaboliteToRemove);

  // Sort all non-input pathways into output bucket and non-participating bucket.
  // NOTE: We only use one thread, so we dedicate one warp to it.
   sortOutputPathways << < 1, 32 >> > (d_binaryVectors, d_metaboliteCoefficients, 
				       pathwayStartIndex + numInputs, // Start sorting after inputs
				       pathwayCount - numInputs,  // Number of pathways to sort
				       metaboliteCount, metaboliteToRemove);
}

void dependencyCheck(int numInputs, int numOutputs, int batch_number) {
   int numBlocks = (numInputs / MAX_THREADS_PER_BLOCK) + 1;
   //printf("dependencyCheck: numInputs=%d, numOutputs=%d, batch_number=%d, numBlocks=%d, batchSize=%d\n", numInputs, numOutputs, batch_number, numBlocks, batchSize);
   dependencyCheck << < numBlocks, MAX_THREADS_PER_BLOCK >> > (d_binaryVectors, d_combinationBins, batchSize, 
           numInputs, //number of inputs
           numInputs + (batch_number * batchSize), //start of next batch of outputs
           numInputs + numOutputs, //start of non-participating
           pathwayCount, pathwayStartIndex);
}

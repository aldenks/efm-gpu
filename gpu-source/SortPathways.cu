#include "SortPathways.h"

__global__
void sortInputPathways(BinaryVector *reactions, float *metaboliteCoefficients,
		  int startIndex, int numToSort, int numReactions, 
		  int numberOfMetabolites, int metaboliteToRemove){

  int tid = blockIdx.x + blockDim.x + threadIdx.x;

  if ( tid == 0 ) {
    int start = startIndex;
    int end = startIndex + numToSort;

    //While the pointers do not overlap
    while(start < end){
      bool is_input_1 = metaboliteCoefficients[start * numberOfMetabolites + 
					       metaboliteToRemove] < NEG_ZERO;
      bool is_input_2 = metaboliteCoefficients[end * numberOfMetabolites + 
					       metaboliteToRemove] < NEG_ZERO;
      if(is_input_1) {
	//Skip this one
	start++;
      } else if(is_input_2){
	//swap the two reactions
	BinaryVector temp = reactions[end];
	reactions[end] = reactions[start];
	reactions[start] = temp;
	
	//swap the two metaboliteCoefficients for pathways
	float tempCoefficient;
	for(int i = 0; i < numberOfMetabolites; ++i){
	  tempCoefficient = metaboliteCoefficients[end * numberOfMetabolites + i];
	  metaboliteCoefficients[end * numberOfMetabolites + i] = 
	    metaboliteCoefficients[start * numberOfMetabolites + i];
	  metaboliteCoefficients[start * numberOfMetabolites + i] = tempCoefficient;
	}
	//move forward
	start++;
	end++;
      } else {
	//Not an input, don't care
	end++;
      }
    }
  }
}

__global__
void sortOutputPathways(BinaryVector *reactions, float *metaboliteCoefficients,
		  int startIndex, int numToSort, int numReactions, 
		  int numberOfMetabolites, int metaboliteToRemove){

  int tid = blockIdx.x + blockDim.x + threadIdx.x;

  if ( tid == 0 ) {
    int start = startIndex;
    int end = startIndex + numToSort;

    //While the pointers do not overlap
    while(start < end){
      bool is_output_1 = metaboliteCoefficients[start * numberOfMetabolites + 
					       metaboliteToRemove] > ZERO;
      bool is_output_2 = metaboliteCoefficients[end * numberOfMetabolites + 
					       metaboliteToRemove] > ZERO;
      if(is_output_1) {
	//Skip this one
	start++;
      }	if(is_output_2){
	//swap the two reactions
	BinaryVector temp = reactions[end];
	reactions[end] = reactions[start];
	reactions[start] = temp;
	
	//swap the two metaboliteCoefficients for pathways
	float tempCoefficient;
	for(int i = 0; i < numberOfMetabolites; ++i){
	  tempCoefficient = metaboliteCoefficients[end * numberOfMetabolites + i];
	  metaboliteCoefficients[end * numberOfMetabolites + i] = 
	    metaboliteCoefficients[start * numberOfMetabolites + i];
	  metaboliteCoefficients[start * numberOfMetabolites + i] = tempCoefficient;
	}
	//move forward
	start++;
	end++;
      } else {
	//Not an output, don't care
	end++;
      }
    }
  }
}

__global__ 
void dependencyCheck(BinaryVector *reactions, int *bins, int batch_size, 
		     int num_inputs, int output_start, int non_part_start,
		     int pathwayCounts){
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int count = 0;
  if(tid < num_inputs){
    BinaryVector input = reactions[tid];
    BinaryVector output, combo, pathway;
    bool is_unique_and_independent = false;
    for(int i = 0; i < batch_size; ++i){
      output = reactions[output_start + i];
      combo = input | output;
      for (int j = 0; is_unique_and_independent && j < pathwayCounts; ++j){
	if(j == tid) continue; //skip this input
	if(j == output_start + i) continue; //skip this output
	pathway = reactions[j];
	if(pathway == combo){
	  //TODO: how can we prevent duplicates?
	  is_unique_and_independent = false;
	  break;
	}
	is_unique_and_independent = ((combo & pathway) != pathway);
      }

      if(is_unique_and_independent){
	bins[batch_size * (i+1) + tid] = i;
	count++;
      } else {
	bins[batch_size * (i+1) + tid] = -1;
      }
    }

    //Assumes first item is 0 to start
    if(output_start == num_inputs)
      bins[tid] = 0;

    bins[tid] += count;
  }
}

void sortInputsOutputs(float *d_metaboliteCoefficients, int pathwayCounts, 
		       BinaryVector *d_reactions, int metaboliteCount, int numInputs, 
		       int numOutputs, int metaboliteToRemove){
  //call the kernel on inputs and outputs
  int numBlocks = (pathwayCounts / MAX_THREADS_PER_BLOCK ) + 1;
  sortInputPathways <<< numBlocks, MAX_THREADS_PER_BLOCK >>> 
    (d_reactions, d_metaboliteCoefficients, 
     0, pathwayCounts, pathwayCounts,
     metaboliteCount, metaboliteToRemove);

  numBlocks = ((pathwayCounts - numInputs) / MAX_THREADS_PER_BLOCK ) + 1;
  sortOutputPathways <<< numBlocks, MAX_THREADS_PER_BLOCK >>> 
    (d_reactions, d_metaboliteCoefficients, 
     numInputs, pathwayCounts-numInputs, pathwayCounts,
     metaboliteCount, metaboliteToRemove);

}

void dependencyCheck(int numInputs, int numOutputs, int batch_number){
  int numBlocks = (numInputs / MAX_THREADS_PER_BLOCK ) + 1;
  dependencyCheck <<< numBlocks , MAX_THREADS_PER_BLOCK >>> 
    (d_binaryVectors, d_combinationBins, batchSize, numInputs,
     numInputs + (batch_number * batchSize), //start of next batch of outputs
     numOutputs, //start of non-participating
     pathwayCount);
}

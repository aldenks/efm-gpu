
__global__
void sortPathways(int *reactions, float *metaboliteCoefficients, int sorting,
		  int startIndex, int numToSort, int numReactions, 
		  int numberOfMetabolites, int metaboliteToRemove){

  int tid = blockIdx.x + blockDim.x + threadIdx.x;

  if(tid < numToSort){
    // Step 1) Each thread searches for an input
    
    //Hope these two go into registers
    int counter = 0;
    int index_of_input = -1;

    //for each pathway
    for(int i = startIndex; i < numReactions; ++i){
      //NOTE: how can we generalize sort w/o divergence
      //REMEMBER: logical-and is slow due to short-circuiting!
      //IDEA: just make 2 kernels? one input, other output

      // if found what I am looking for, counter++;
      // NOTE: Inputs are determined by negative coefficents, 
      //       outputs are positive

      if(sorting < 0) { //Search for inputs
	if(metaboliteCoefficients[i * numberOfMetabolites + 
				  metaboliteToRemove] < 0) {
	  counter++;
	}
      } else { //Search for outputs (Assuming sorting != 0)
	if(metaboliteCoefficients[i * numberOfMetabolites + 
				  metaboliteToRemove] > 0) {
	  counter++;
	}
      }
      
      //if (counter - 1) = (blockIdx.x * blockDim.x + tid)
      if( (counter-1) == tid ){
	//  then index_of_input = current_index      
	index_of_input = i;
      }
    }
  // end for

  // syncthreads to prevent writing while other threads are reading
    __syncthreads();

  // swap (blockIdx.x * blockDim.x + tid) slot with index_of_input
  // NOTE: will need to copy data to LOCAL memory
  // QUESTION: how do we know we will NOT have a write/read conflict 
  //           between threads?
  //           i.e.: thread 0 is swapping with index 3 while
  //                 thread 3 is swapping with index 7
  // also swap metabolite coefficients
  // NOTE: each row is one pathway
  }
}

void sortInputsOutputs(float *d_metaboliteCoefficients, int pathwayCounts, 
		       int *reactions, int metaboliteCount, int numInputs, 
		       int numOutputs, int metaboliteToRemove){
  //call the kernel on inputs and outputs
}

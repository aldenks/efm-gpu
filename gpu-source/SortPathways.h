#ifndef SORT_PATHWAYS_H
#define SORT_PATHWAYS_H

#include "Setup.h"

/*
 * Takes the circular buffer of reactions and metabolites d_reactions
 * and sorts them into buckets of inputs, outputs, and non-participating
 * according to a chosen metabolite
 *
 * +------------------------------------------------------------------+
 * |                                                                  |
 * +------------------------------------------------------------------+
 *
 * +------------------------------------------------------------------+
 * |    inputs        |     outputs        |     non-participating    |
 * +------------------------------------------------------------------+
 *
 * Also takes the metabolite coefficients 2D array and sorts each row
 * into the same buckets. Asumming we have m metabolites and n
 * pathways, we have
 *
 *       m1   m2   m3   m4   ...  mr  ...   mm
 *     +----+----+----+----+----+----+----+----+
 * r1  |    |    |    |    |    | 0  |    |    |
 *     +----+----+----+----+----+----+----+----+
 * r2  |    |    |    |    |    | 2  |    |    |
 *     +----+----+----+----+----+----+----+----+
 * r3  |    |    |    |    |    | -1 |    |    |
 *     +----+----+----+----+----+----+----+----+
 * r4  |    |    |    |    |    | -3 |    |    |
 *     +----+----+----+----+----+----+----+----+
 * ... |    |    |    |    |    |    |    |    |
 *     +----+----+----+----+----+----+----+----+
 * rn  |    |    |    |    |    | 0  |    |    |
 *     +----+----+----+----+----+----+----+----+
 *
 * become this after sorting
 *
 *       m1   m2   m3   m4   ...  mr  ...   mm
 *     +----+----+----+----+----+----+----+----+  ---+
 * r3  |    |    |    |    |    | -1 |    |    |     |
 *     +----+----+----+----+----+----+----+----+     +--> inputs
 * r4  |    |    |    |    |    | -3 |    |    |     |
 *     +----+----+----+----+----+----+----+----+  ---+
 * r2  |    |    |    |    |    | 2  |    |    |  ------> outputs
 *     +----+----+----+----+----+----+----+----+  ---+
 * r1  |    |    |    |    |    | 0  |    |    |     |
 *     +----+----+----+----+----+----+----+----+     |
 * ... |    |    |    |    |    |    |    |    |     +---> non-part
 *     +----+----+----+----+----+----+----+----+     |
 * rn  |    |    |    |    |    | 0  |    |    |     |
 *     +----+----+----+----+----+----+----+----+  ---+
 */
void sortInputsOutputs(int numInputs, int numOutputs, int metaboliteToRemove);

/*
 * Takes the array of binary vectors ordered into inputs, outputs,
 * and non-participating. Generates combinations by combining one input
 * from inputs with one output from outputs and checks dependency against
 * all binary vectors in this array.
 * +------------------------------------------------------------------+
 * |    inputs        |     outputs        |     non-participating    |
 * +------------------------------------------------------------------+
 *
 * Output is written in bins, where each column is a bin for a particular
 * input. The first row has the number of output indices listed in the bin.
 * The rest of the rows contain output indices.
 *
 *       0    1    2    3    4    5    6   ...  n_in
 *     +----+----+----+----+----+----+----+----+----+
 *     | 3  | 1  | 0  | 0  | 2  | 0  | 0  |    | 2  | ------> counts
 *     +----+----+----+----+----+----+----+----+----+ --+
 *     | 15 | 12 |    |    | 13 |    |    |    | 21 |   |  
 *     +----+----+----+----+----+----+----+----+----+   |
 *     | 17 |    |    |    | 15 |    |    |    | 26 |   |
 *     +----+----+----+----+----+----+----+----+----+   |
 *     | 26 |    |    |    |    |    |    |    |    |   +--> output_indices
 *     +----+----+----+----+----+----+----+----+----+   |
 *     |    |    |    |    |    |    |    |    |    |   |
 *     +----+----+----+----+----+----+----+----+----+   |
 *     |    |    |    |    |    |    |    |    |    |   |
 *     +----+----+----+----+----+----+----+----+----+ --+
 *
 * ASSUMPTIONS: The bin is considered zeroed out at the beginning of
 *              every call to this function.
 *              The number of columns is exactly numInputs.
 *              The number of rows is exactly batchSize.
 */
void dependencyCheck(int numInputs, int numOutputs, int batch_number);


//Debugging functions
void checkSorting(int *h_sorts, int *d_sorts, int metaboliteToRemove);

#endif

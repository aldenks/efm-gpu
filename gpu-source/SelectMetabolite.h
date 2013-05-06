#ifndef SELECT_METABOLITE_H
#define SELECT_METABOLITE_H

int getNextMetabolite(float* d_metaboliteCoefficients, int pathwayStartIndex, int pathwayCount, int metaboliteCount, bool* d_balancedMetabolites, int* d_inputCounts, int* d_outputCounts, int* h_inputCounts, int* h_outputCounts);

#endif	/* SELECT_METABOLITE_H */

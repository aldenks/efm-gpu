#ifndef SELECT_METABOLITE_H
#define	SELECT_METABOLITE_H

int getNextMetabolite(float* d_metaboliteCoefficients, int pathwayCount, int metaboliteCount, bool* d_balancedMetabolites, int* d_inputCounts, int* d_outputCounts, int* h_inputCounts, int* h_outputCounts);

__global__
void computeMetaboliteInputOutputCounts(float* metaboliteCoefficients, int pathwayCount, int metaboliteCount, bool* balancedMetabolites, int* inputCounts, int* outputCounts);

#endif	/* SELECT_METABOLITE_H */

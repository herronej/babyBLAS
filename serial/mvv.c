#ifdef __cplusplus
extern "C" {
#endif
	void mvv_(int *threads, int *N, double* mat, double* vec, double* vresults);
#ifdef __cplusplus
	}
#endif

#include <stdio.h>

void mvv_(int *threads, int *N, double* mat, double* vec, double* vresults){
	int i, j;

	int length = *N;

	for(i=0;i<length;i++){
		for(j=0;j<length;j++){
			*(vresults+i) += *(mat+(length*i)+j) * *(vec+j);
			//double r = *(vresults+i);
			//printf("%f\n", *(vresults+i));


		}
	}
}

#ifdef __cplusplus
extern "C" {
#endif
	void dot_(int *threads, int *N, double *vec1, double *vec2, double *rval);
#ifdef __cplusplus
	}
#endif

#include <stdio.h>

void dot_(int *threads, int *N, double *vec1, double *vec2, double *rval){
	int i;

	int len = *N;

	//double rval = 0.0;

	for(i=0; i < len; i++){

		//printf("%f\n", *(vec1+i) * *(vec2+i));

		*rval += (*(vec1+i) * *(vec2+i));
	}	

	//return rval;
}

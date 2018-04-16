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

	for(i=0; i < len; i++){

		*rval += (*(vec1+i) * *(vec2+i));
	}
}

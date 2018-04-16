#ifdef __cplusplus
extern "C" {
#endif
	void dot_(int *threads, int *N, double *vec1, double *vec2, double *rval);
#ifdef __cplusplus
	}
#endif

#include <omp.h>

void dot_(int *threads, int *N, double *vec1, double *vec2, double *rval){
	int i;

	int len = *N;

    omp_set_num_threads(*threads);
	
	#pragma omp parallel for
	for(i=0; i < len; i++){

                *rval += (*(vec1+i) * *(vec2+i));

		//double x = *rval;

		//printf("rval: %d\n", x);
        }	
}


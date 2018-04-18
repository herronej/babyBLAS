#ifdef __cplusplus
extern "C" {
#endif
	void mvv_(int *threads, int *N, double* mat, double* vec, double* vresults);
#ifdef __cplusplus
	}
#endif

#include <omp.h>

void mvv_(int *threads, int *N, double* mat, double* vec, double* vresults){
	int i, j;

	int length = *N;

	omp_set_num_threads(*threads);

//#pragma omp parallel shared(length) private(i,j)
//{
	/*for(i=0;i<length;i++){
		for(j=0;j<length;j++){
			*(vresults+i) += *(mat+(length*i)+j) * *(vresults+j);
		}
	}*/

	#pragma omp parallel for
	for(i=0;i<length;i++){
                for(j=0;j<length;j++){
                        *(vresults+i) += *(mat+(length*i)+j) * *(vec+j);
                        //printf("%f\n", *(vresults+i));
                }
        }

//}
}

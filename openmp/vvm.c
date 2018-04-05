#ifdef __cplusplus
extern "C" {
#endif
	void vvm_(int *threads, int *len, double *va, double *vb, double *ma);
#ifdef __cplusplus
	}
#endif

#include <omp.h>

void vvm_(int *threads, int *len, double *va, double *vb, double *ma){
	int i, j;
	int alength = *len;

	omp_set_num_threads(*threads);

#pragma omp parallel shared(alength) private(i,j)
{

	for(i=0;i<alength;i++){
		for(j=0;j<alength;j++){
			*(ma+(alength*i)+j) = *(va+i) * *(vb+j);
		}
	}
}
}


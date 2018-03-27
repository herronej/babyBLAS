#ifdef __cplusplus
extern "C" {
#endif
	void mvv_(int *N, double* mat, double* vec, double* vresults);
#ifdef __cplusplus
	}
#endif

void mvv_(int *N, double* mat, double* vec, double* vresults){
	int i, j;

	int length = *N;

	for(i=0;i<length;i++){
		for(j=0;j<length;j++){
			*(vresults+i) += *(mat+(length*i)+j) * *(vresults+j);
		}
	}
}

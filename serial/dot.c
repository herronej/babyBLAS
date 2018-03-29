#ifdef __cplusplus
extern "C" {
#endif
	void dot_(int *N, double *vec1, double *vec2, double *rvec);
#ifdef __cplusplus
	}
#endif

void dot_(int *N, double *vec1, double *vec2, double *rvec){
	int i;

	int len = *N;

	for(i=0; i < len; i++){
		*(rvec+i) = *(vec1+i) + *(vec2+i);
	}	
}

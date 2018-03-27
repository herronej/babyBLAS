#ifdef __cplusplus
extern "C" {
#endif
	void dot_(int *N, double *vec1, double *vec2);
#ifdef __cplusplus
	}
#endif

void dot_(int *N, double *vec1, double *vec2, double *rvec){
	int i;
	for(i=0; i < N; i++){
		*(rvec+i) = *(rvec+i) + *(rvec+j);
	}	
}

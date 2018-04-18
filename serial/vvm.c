#ifdef __cplusplus
extern "C" {
#endif
	void vvm_(int *threads, int *len, double *va, double *vb, double *ma);
#ifdef __cplusplus
	}
#endif

void vvm_(int *threads, int *len, double *va, double *vb, double *ma){
	int i, j;
	int alength = *len;

	#ifdef STRIP8
        const int stride = 8;
        int mod = alength % stride;

        for(i=0;i<alength;i++){

		
            	for(j=0;j<mod;j++){
                        *(ma+(alength*i)+j) = *(va+i) * *(vb+j);
                }
		for(j=0;j<alength;j+=stride){
                        *(ma+(alength*i)+(j)) = *(va+i) * *(vb+j);
			*(ma+(alength*i)+(j+1)) = *(va+i) * *(vb+j+1);
			*(ma+(alength*i)+(j+2)) = *(va+i) * *(vb+j+2);
			*(ma+(alength*i)+(j+3)) = *(va+i) * *(vb+j+3);
			*(ma+(alength*i)+(j+4)) = *(va+i) * *(vb+j+4);
			*(ma+(alength*i)+(j+5)) = *(va+i) * *(vb+j+5);
			*(ma+(alength*i)+(j+6)) = *(va+i) * *(vb+j+6);
			*(ma+(alength*i)+(j+7)) = *(va+i) * *(vb+j+7);
                }
        }
	#else
	for(i=0;i<alength;i++){
		for(j=0;j<alength;j++){
			*(ma+(alength*i)+j) = *(va+i) * *(vb+j);
		}
	}
	#endif
}

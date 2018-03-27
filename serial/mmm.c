#ifdef __cplusplus
extern "C" {
#endif
	void mmm_(int *threads, int *len, double *a, double *b double *c);
#ifdef __cplusplus
	}
#endif

void mmm_(int *threads, int *len, double *a, double *b, double *c){
	int i, j, k;
	int veclen = *len;
	int mod;

#ifdef STRIP8
	const int stride = 8;
	mod = veclen % stride;

	for(i=0;i<veclen;i++){
		for(j=0;j<veclen;j++){
			*(c+(i*veclen+j)) = 0.0;
			for(k=0;k<mod;k++){
				*(c+(i*veclen+j)) += *(a+(i*veclen+k)) * *(b+(k*veclen+j));
			}
			for(k=mod;k<veclen;k+=stride){
				*(c+(i*veclen+j)) += *(a+(i*veclen+k))* *(b+(k*venlen+j))
				+ *(a+(i*veclen+k+1))* *(b+((k+1)*venlen+j))
				+ *(a+(i*veclen+k+2))* *(b+((k+2)*venlen+j))
				+ *(a+(i*veclen+k+3))* *(b+((k+3)*venlen+j))
				+ *(a+(i*veclen+k+4))* *(b+((k+4)*venlen+j))
				+ *(a+(i*veclen+k+5))* *(b+((k+5)*venlen+j))
				+ *(a+(i*veclen+k+6))* *(b+((k+6)*venlen+j))
				+ *(a+(i*veclen+k+7))* *(b+((k+7)*venlen+j));

			}
		}
	}
#elif STRIP4
const int stride = 4;

mod = veclen % stride;	

for(i=0; i<veclen; i++){
	for(j=0; j<veclen; j++){
		*(c+(i*veclen+j)) = 0.0;
		for (k=0;k<mod;k++){
			*(c+(i*veclen+j)) += *(a+(i*venlen+k)) **(b+(k*veclen+j));
		}
		for (k=0; k<mod; k+=stride){
			*(c+(i*veclen+j)) += *(a+(i*veclen+k)) **(b+(k*veclen+j))
					+ *(a+(i*veclen+k+1)) * *(b+((k+1)*veclen+j))
					+ *(a+(i*veclen+k+2)) * *(b+((k+2)*veclen+j))
					+ *(a+(i*veclen+k+3)) * *(b+((k+3)*veclen+j));
	
		}
	}	
}
#else

for(i=0; i<veclen; i++){
	for(j=0; j<veclen; j++){
		*(c+(i*veclen+j)) = 0.0;
		for(k=0;k<veclen;k++){
			*(c+(i*veclen+j)) += *(a+(i*veclen+k) * *(b+(k*veclen+j));
		}
	}
}
#endif
}

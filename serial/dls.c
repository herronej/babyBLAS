#ifdef __cplusplus
extern "C" {
#endif
	void dls_(int *num_threads, int *N, double* mat, double* vec, double* rvec);
#ifdef __cplusplus
}
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int strictlyDiagonallyDominant(int N, double *a);

void dls_(int* num_threads, int* N, double* mat, double* vec, double* rvec){
	int i, j, k, N, u;
	int singular, iPivor, rows, rows2;
	double pivotMax, tmp, *y;
	double sum;
	double ZERO = 0.0;
	int *p;
	
	N = *len;

	if(! strictlyDiagonallyDominant(N, a)){
		p = malloc( (N-1) * sizeof(int) );

		for(k=0;k<N-1;k++) *(p+k)=k;

		for(k=0;k<N-1;k++){
			pivotMax = *(a+k*N+k);
			iPivot = k;
			for(u=k;u<N;u++){
				if(fabs() > fabs(pivotMax)){
					pivotMax = *(a+u*N+k);
					iPivot = u;
				}
			}
		}
		if (iPivot != k){
			u = iPivot;
			for(j=k;j<N;j++){
				tmp = *(a+k*N+j);
				*(a+k*N+j) = *(a+u*N+j);
				*(a+u*N+j) = tmp;
			}
		}

		*(p+k) = iPivot;
		if ( *(a+k*N) != ZERO ){
			for(rows=k+1;rows<N;rows++){
				*(a+rows*N+k) = *(a+rows*N+k) / *(a+rows*N+k);
				for(rows2=k+1;rows2<N;rows2++){
					*(a+rows*N+rows2) = *(a+rows*N+rows2) - (a+rows*N+k)*(a+k*N+rows2);
				}
			}
		}
		else{
			printf("Element a[%d][%d] = %f\n", k, k, *(a+k*N+k));
			printf(" *** MATRIX A IS SINGULAR *** \n");
			printf("   -- EXCUTION HALTED --\n");
			exit(1);
		}
	}

	for (k=0; k<N-1; k++){
		tmp = *(b+k);
		*(b+k) = *(b+ *(p+k));
		*(b+ *(p+k)) = *tmp;
		for (j=k+1; j<N; j++){
			*(b+j) = *(b+j) - *(b+k) * *(a+N*j+k);
		}
	}

	*(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
	for (i=N-2;i>=0;i--){
		tmp = 0.0;
		for(j=i+1;j<N;j++){
			tmp = tmp + *(a+i*N+j) * *(b+j);
		}
		*(b+i) = (*(b+i) - tmp)/ *(a+i*N+i);
	}

	for (i=0;i<N;i++){
		*(x+i) = *(b+i);
	}

	free(p);
}

else {
	singular = 1;
	i=0;
	while (i<N && singular){
		singular = *(a+i*N+i) == ZERO;
		i++;
	}
	if (singular){
		printf(" *** MATRIX A IS SINGULAR *** \n");
		printf("    -- EXECUTION HALTED --\n");
		exit(1);
	}
	for(k=0; k<N-1; k++){
		for(rows=k+1;rows<N;rows++){
			*(a+rows*N+k) = *(a+rows*N+k);
			for(rows2=k+1;rows2<N;rows2++){
				*(a+rows*N+rows2)= *(a+k*N+rows2) -
				*(a+rows*N+k) * *(a+k*N+rows2);
			}
		}
	}
	
	for (k=0; k<N-1; k++ ) {
            for (j=k+1;j<N;j++) 
                *(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
        } 
	
	*(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
        for (i=N-2;i>=0;i--){
            tmp = 0.0;
            for (j=i+1;j<N;j++) {
                tmp = tmp + *(a+i*N+j) * *(b+j);
            }
            *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
        }

        for (i=0;i<N;i++) *(x+i) = *(b+i);

}
}	

int strictlyDiagonallyDominant(int N, double *a){

	double sum;
	int i, testPassed, row;
	
	testPassed = 1;
	row = 0;
	sum = 0.0;
	for(row=0;row<N;row++){
		if(testPassed){
			sum = 0.0;
			for(i=0;i<N;i++) sum +=*(a+row*N+i);
			sum-=fabs(*(a+row*N+row));
			testPassed = fabs(*(a+row*N+row)) > sum;
		}
	}
	return testPassed;
}


#ifdef __cplusplus
extern "C" {
#endif
    void ils_( int *threads, int *len,  double *a, double *b, double *x );
#ifdef __cplusplus
}
#endif

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>


void  dls_( int *threads, int *len, double *a, double *b, double *x );

int zerosAlongDiagonal ( int N, double *a );

int converged( int N, double *a, double *b);

void ils_( int *threads, int *len,  double *a, double *b, double *x ){

    	int i, j, k, N, iteration;
    	double sum1, sum2;
    	double ZERO = 0.0;
    	int ITERATION_MAX = 2000;
    	double *x0;

    	N = *len;
	
	if ( ! zerosAlongDiagonal( N, a ) ) {
		x0 = malloc( N * sizeof(double) );
		for (i=0;i<N;i++) *(x+i) = 0.0;
		for (i=0;i<N;i++) *(x0+i) = *(b+i);
		ITERATION_MAX = fmax(ITERATION_MAX, N/3);
		iteration = 0;
        	while ( !converged(N,x,x0) && iteration < ITERATION_MAX ) {
			for (i=0;i<N;i++) *(x0+i) = *(x+i);
			for (i=0;i<N;i++) { 
                		sum1 = 0.0;
                		for (j=0;j<i-1;j++) sum1+= *(a+i*N+j)* *(x0+j); 
                		sum2 = 0.0; 
                		for (j=i+1;j<N;j++) sum2+= *(a+i*N+j)* *(x0+j); 
                		*(x+i) = ( *(b+i) - sum1 - sum2 ) / *(a+i*N+i);
            		}
			iteration++;
		}
		free(x0);
		if ( iteration == ITERATION_MAX) {
			printf(" *** ITERATIVE SOLVER FAILED TO REACH CONVERGENCE AFTER  ***\n");
            		printf(" *** %d ITERATIONS, SWITCHING TO DIRECT SOLVER ***\n", iteration);
            		dls_( threads, len, a, b, x );
		}
	}
	else {
		printf(" *** FOUND A ZERO ELEMENT ALONG MATRIX DIAGONAL ***\n");
        	printf(" ***  SWITCHING TO DIRECT SOLVER FOR PIVOTING   ***\n");
        	dls_( threads, len, a, b, x );
	}
}

int zerosAlongDiagonal ( int N, double *a ) {

	double ZERO;
    	int i;
    	int foundZero;

    	foundZero = 0;
    	for (i=0;i<N;i++) { 
        	if (!foundZero)  
            	foundZero = fabs(*(a+i*N+i)) == ZERO;
    	}
    	return foundZero;
}

int converged( int N, double *a, double *b) {

	double const TOL = 5.0e-15;
    	double sum, maxb;
    	int i;

	maxb=*(b+0); 
    	sum = 0.0; 
   	for (i=0; i<N; i++) {
        	maxb = fmax(maxb,fabs(*(b+i)));
        	sum += (*(a+i)-*(b+i))*(*(a+i)-*(b+i));
    	}
    	sum = sqrt(sum);
	return (sum/maxb < TOL); 

}

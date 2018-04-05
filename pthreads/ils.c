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
#include <pthread.h>

void  ils_( int *threads, int *len, double *a, double *b, double *x );

int zerosAlongDiagonal ( int N, double *a );

int converged( int N, double *a, double *b);

void *ils_thread_worker();

struct args {
	int N;
	int n_threads;
	int startRow;
	int stopRow;
	double *aptr;
	double *bptr;
	double *xptr;
};

void ils_( int *threads, int *len,  double *a, double *b, double *x ){

	int numThreads = *threads;
    	int matrixDimension = *len;
    	int *numberOfRows;
    	int startRow, stopRow;
    	pthread_t *thread_id;
    	struct args *thread_args;

    	int i, j, k, N, iteration;
    	double sum1, sum2;
    	double ZERO = 0.0;
    	int ITERATION_MAX = 2000;
    	double *x0;

    	N = *len;
	
	if( matrixDimension < numThreads){
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

	else{

		thread_id = (pthread_t *) malloc(numThreads * sizeof(pthread_t));
		numberOfRows = (int*) malloc(numThreads * sizeof(int));

		for(int i = 0; i < numThreads; i++){
			*(numberOfRows+i) = matrixDimension / numThreads;
		}

		for(int i = 0; i< matrixDimension % numThreads; i++){
			*(numberOfRows+i) = *(numberOfRows+i) + 1;
		}

		stopRow = 0;
		for(int i = 0; i < numThreads; i++){
			{
				startRow=stopRow;
                		stopRow=startRow+*(numberOfRows+i);
                		thread_args = ( struct args * )  malloc(sizeof( struct args));
                		thread_args->N   = matrixDimension;
				thread_args->n_threads = threads;
                		thread_args->startRow = startRow;
                		thread_args->stopRow = stopRow; 
                		thread_args->aptr = a;
                		thread_args->bptr = b;
                		thread_args->xptr = x;

                		pthread_create( thread_id+i, NULL, &ils_thread_worker, thread_args );
		
			}
	
		}

		for(int i = 0; i < numThreads; i++){
			pthread_join( *(thread_id+i), NULL);
		}

		free(numberOfRows);
		free(thread_id);
	}

}

void *ils_thread_worker(struct args *thread_args){

	int i, j, k;
	double val;	
	int rowStart, rowStop, N;
	double *a, *b, *x;

	N = thread_args->N;
    	rowStart = thread_args->startRow;
    	rowStop = thread_args->stopRow; 
	int *threads = thread_args->n_threads;
    	a = thread_args->aptr;
    	b = thread_args->bptr;
    	x = thread_args->xptr;

	int iteration;
    	double sum1, sum2;
    	double ZERO = 0.0;
    	int ITERATION_MAX = 2000;
    	double *x0;

	if ( ! zerosAlongDiagonal( N, a ) ) {
		x0 = malloc( N * sizeof(double) );
		for (i=rowStart;i<rowStop;i++) *(x+i) = 0.0;
		for (i=rowStart;i<rowStop;i++) *(x0+i) = *(b+i);
		ITERATION_MAX = fmax(ITERATION_MAX, N/3);
		iteration = 0;
        	while ( !converged(N,x,x0) && iteration < ITERATION_MAX ) {
			for (i=rowStart;i<rowStop;i++) *(x0+i) = *(x+i);
			for (i=rowStart;i<rowStop;i++) { 
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
            		dls_( threads, N, a, b, x );
		}
	}
	else {
		printf(" *** FOUND A ZERO ELEMENT ALONG MATRIX DIAGONAL ***\n");
        	printf(" ***  SWITCHING TO DIRECT SOLVER FOR PIVOTING   ***\n");
        	dls_( threads, N, a, b, x );
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

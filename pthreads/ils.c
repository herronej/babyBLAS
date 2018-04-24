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

struct args6 {
	int N;
        int startRow;
        int stopRow;
        double *Aptr;
        double* X0ptr;
        double* Xptr;
        double* Bptr;
};

struct args7 {
        int N;
        int startRow;
        int stopRow;
        double* X0ptr;
        double* Xptr;
        double* Bptr;
};


void *dls_thread_worker7( struct args6 *thread_args  );
void *dls_thread_worker8( struct args7 *thread_args  );


// Need the following prototype in case of a zero along the diagonal
void  dls_( int *threads, int *len, double *a, double *b, double *x );

// Prototype for code to check for zeros along the diagonal
int zerosAlongDiagonal ( int N, double *a );

// Prototype for code to check for convergence
int converged( int N, double *a, double *b);

void ils_( int *threads, int *len,  double *a, double *b, double *x ){

    /* in serial code, *threads not used. It is retained here so the code can be called
     * identically to the threaded methods.
     */


    int i, j, k, N, iteration;
    double sum1, sum2;
    double ZERO = 0.0;
    int ITERATION_MAX = 2000;
    double *x0;

    int numThreads = *threads;
    int *numberOfRows;
    int startRow, stopRow;
    pthread_t *thread_id;
    struct args6 *thread_args6;
    struct args7 *thread_args7;

    N = *len;

    // The iterative Jacobi method will fail if there are zeros along the diagonal of the
    // matrix.  We could reorder the matrix, but in this case if zeros are found we will switch
    // to our discrete linear solver which will do a form of Gaussian Elimination with partial
    // pivoting to solve the system. 

    // NOTE: we are not modifying matrix A or vector B in any manner in this iterative procedure
    // so that if we fail to achieve convergence we can drop back to the direct solver with the original
    // A matrix and B vector.
    if ( ! zerosAlongDiagonal( N, a ) ) {

        // Do Jacobi Iterative Method to solve Ax=b. 

        // Create a temporary vector to hold initial values and intermediate steps 

        x0 = malloc( N * sizeof(double) );

        // Fill the x0 vector with initial values of zero
	if(N < numThreads){
	        for (i=0;i<N;i++) *(x+i) = 0.0;

        	for (i=0;i<N;i++) *(x0+i) = *(b+i);
        }
	else{
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
                numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

                for (int i=0; i<numThreads; i++ ){
                        *(numberOfRows+i) = (N-k-1) / numThreads;
                }
                for (int i=0; i< (N-k-1) % numThreads; i++ ){
                        *(numberOfRows+i) = *(numberOfRows+i) + 1;
                }

                stopRow=0;
                for(int i=0; i < numThreads ; i++) {
                   {
                        startRow=stopRow;
                        stopRow=startRow+*(numberOfRows+i);
                        thread_args7 = ( struct args7 * )  malloc(sizeof( struct args7));
                        thread_args7->N   = N;
                        thread_args7->startRow = startRow;
                        thread_args7->stopRow = stopRow;
                        thread_args7->X0ptr = x0;
                        thread_args7->Xptr = x;
                        thread_args7->Bptr = b;

                        pthread_create( thread_id+i, NULL, &dls_thread_worker8, thread_args7 );
                   }
                }
                for(int i=0; i < numThreads ; i++) {
                        pthread_join( *(thread_id+i), NULL);
                }

                free(numberOfRows);
                free(thread_id);
	}
        // If more than N/3 iterations are done, the direct solver is more efficient
        ITERATION_MAX = fmax(ITERATION_MAX, N/3);

        iteration = 0;
        while ( !converged(N,x,x0) && iteration < ITERATION_MAX ) {

            // copy last result to initial values

            for (i=0;i<N;i++) *(x0+i) = *(x+i);

            // start the reduction process  (ref: Golub and van Loan, Chapter 10)
            if(N < numThreads){
	            for (i=0;i<N;i++) { 
        	        sum1 = 0.0;
                	for (j=0;j<i-1;j++) sum1+= *(a+i*N+j)* *(x0+j); 
                	sum2 = 0.0; 
                	for (j=i+1;j<N;j++) sum2+= *(a+i*N+j)* *(x0+j); 
                	*(x+i) = ( *(b+i) - sum1 - sum2 ) / *(a+i*N+i);
            	    }
	    }
	    else{
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
		numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

		for (int i=0; i<numThreads; i++ ){
            		*(numberOfRows+i) = (N-k-1) / numThreads;
       		}
        	for (int i=0; i< (N-k-1) % numThreads; i++ ){
           		*(numberOfRows+i) = *(numberOfRows+i) + 1;
       		}

		stopRow=k+1;
        	for(int i=0; i < numThreads ; i++) {
            	   {   
                	startRow=stopRow;
                	stopRow=startRow+*(numberOfRows+i);
                	thread_args6 = ( struct args6 * )  malloc(sizeof( struct args6));
            		thread_args6->N   = N;
                	thread_args6->startRow = startRow;
                	thread_args6->stopRow = stopRow; 
      	 		thread_args6->Aptr = a;
			thread_args6->X0ptr = x0;
			thread_args6->Xptr = x;
			thread_args6->Bptr = b;

            		pthread_create( thread_id+i, NULL, &dls_thread_worker7, thread_args6 );
            	   }
        	}
        	for(int i=0; i < numThreads ; i++) {
      			pthread_join( *(thread_id+i), NULL); 
        	}

        	free(numberOfRows);
        	free(thread_id);

	    }
            iteration++;

        }

        // the initial value array is no longer needed
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


// Code to check for zeros along the diagonal
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

// Code to check for convergence (A. Pounds, 2018)
int converged( int N, double *a, double *b) {

    // Compute the distance between the vectors and see if the 2-Norm is
    // within tolerance

    double const TOL = 5.0e-15;
    double sum, maxb;
    int i;

    // find max in array b for tolerance scaling while computing sum

    maxb=*(b+0); 
    sum = 0.0; 
    for (i=0; i<N; i++) {
        maxb = fmax(maxb,fabs(*(b+i)));
        sum += (*(a+i)-*(b+i))*(*(a+i)-*(b+i));
    }
    sum = sqrt(sum);

    // by dividing by the largest value in the b matrix we effectively
    // scale the 2-Norm so it can achieve machine precision
    return (sum/maxb < TOL);    

}

void *dls_thread_worker7( struct args6 *thread_args  ) {

    int i, j, k;
    double val, tmp, sum;
    int sum1, sum2, rowStart, rowStop, N;
    double *a, *x0, *x, *b;

    N = thread_args->N;
    rowStart = thread_args->startRow;
    rowStop = thread_args->stopRow;
    a = thread_args->Aptr;
    x0 = thread_args->X0ptr;
    x = thread_args->Xptr;
    b = thread_args->Bptr;


    //printf("rowStart: %d rowStop: %d\n", rowStart, rowStop);

    for (i=rowStart;i<rowStop;i++) {
        sum1 = 0.0;
        for (j=0;j<i-1;j++) sum1+= *(a+i*N+j)* *(x0+j);
        sum2 = 0.0;
        for (j=i+1;j<N;j++) sum2+= *(a+i*N+j)* *(x0+j);
    	*(x+i) = ( *(b+i) - sum1 - sum2 ) / *(a+i*N+i);
    }

    free(thread_args);
    pthread_exit(NULL);
}

void *dls_thread_worker8( struct args7 *thread_args  ) {

    int i, j, k;
    double val, tmp, sum;
    int sum1, sum2, rowStart, rowStop, N;
    double *a, *x0, *x, *b;

    N = thread_args->N;
    rowStart = thread_args->startRow;
    rowStop = thread_args->stopRow;
    x0 = thread_args->X0ptr;
    x = thread_args->Xptr;
    b = thread_args->Bptr;


    //printf("rowStart: %d rowStop: %d\n", rowStart, rowStop);

    for (rowStart=0;i<rowStop;i++) *(x+i) = 0.0;

    for (rowStart=0;i<rowStop;i++) *(x0+i) = *(b+i);


    free(thread_args);
    pthread_exit(NULL);
}

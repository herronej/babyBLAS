#ifdef __cplusplus
extern "C" {
#endif
    void dls_( int *threads, int *len,  double *a, double *b, double*c );
#ifdef __cplusplus
}
#endif

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include <pthread.h>

struct args {
    int N;
    int k;
    int startRow;
    int stopRow;
    double *Aptr;
};

struct args2 {
    int N;
    int startRow;
    int stopRow;
    double *Aptr;
    double *Bptr;
};

struct args3 {
    int N;
    int startRow;
    int stopRow;
    double *Aptr;
    double *Bptr;
    int *Pptr;
};

struct args4 {
    int N;
    int startRow;
    int stopRow;
    int *Pptr;
};

struct args5 {
    int startRow;
    int stopRow;
    double *Xptr;
    double *Bptr;
};


void *dls_thread_worker( struct args *thread_args  );
void *dls_thread_worker2( struct args2 *thread_args  );
void *dls_thread_worker3( struct args3 *thread_args  );
void *dls_thread_worker4( struct args *thread_args  );
void *dls_thread_worker5( struct args4 *thread_args  ); 
void *dls_thread_worker6( struct args5 *thread_args  );

/* Function prototype for code used in dls */
int strictlyDiagonallyDominant( int N, double *a ); 


void dls_( int *threads, int *len,  double *a, double *b, double *x ){

    /* in serial code, *threads not used. It is retained here so the code can be called
     * identically to the threaded methods.
     */


    int i, j, k, N, u;
    int singular, iPivot, rows, rows2;
    double pivotMax, tmp, *y;
    double sum;
    double ZERO = 0.0;
    int *p;

    int numThreads = *threads;
    int *numberOfRows;
    int startRow, stopRow;
    pthread_t *thread_id;
    struct args *thread_args;
    struct args2 *thread_args2;	
    struct args3 *thread_args3;
    struct args4 *thread_args4;
    struct args5 *thread_args5;

    N = *len;

    // Check A for strict diagonal dominance to see if we can reduce the matrix without 
    // doing any row interchanges.   We could also check for positive definiteness to
    // achieve the same thing.

    if ( ! strictlyDiagonallyDominant( N, a ) ) {

        // Do Gaussian Elimination with Partial Pivoting 
        //   (modified from Golub and van Load, Chapter 3) 

        // Create an array to hold pivot swaps 

        p = malloc( (N-1) * sizeof(int) );

	if(N < numThreads){
        	for (k=0;k<N-1;k++) *(p+k)=k;
	}
	else{
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
                        numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

                        for (int i=0; i<numThreads; i++ ){
                	        *(numberOfRows+i) = (N-1) / numThreads;
                        }
                        for (int i=0; i< (N-1) % numThreads; i++ ){
                                *(numberOfRows+i) = *(numberOfRows+i) + 1;
                        }

                        stopRow=0;
                        for(int i=0; i < numThreads ; i++) {
                            {
                                startRow=stopRow;
                                stopRow=startRow+*(numberOfRows+i);
                                thread_args4 = ( struct args4 * )  malloc(sizeof( struct args4));
                                thread_args4->N   = N;
                                thread_args4->startRow = startRow;
                                thread_args4->stopRow = stopRow;
                                thread_args4->Pptr = p;

                         	pthread_create( thread_id+i, NULL, &dls_thread_worker5, thread_args4 );
                             }
                         }
                         for(int i=0; i < numThreads ; i++) {
                         	pthread_join( *(thread_id+i), NULL);
                         }

                         free(numberOfRows);
                         free(thread_id);
	}

        // Search for largest value in the column and swap the 
        // entire row containing that value with the current
        // pivot row.

        for (k=0;k<N-1;k++) {
            pivotMax = *(a+k*N+k);
            iPivot = k; 
            for (u=k;u<N;u++) {
                if ( fabs(*(a+u*N+k)) > fabs(pivotMax) ) {
                    pivotMax = *(a+u*N+k);
                    iPivot = u;
                }
            }
            // If a greater pivot value was found, swap the rows.
            if ( iPivot != k ) {
                u = iPivot; 
                for (j=k;j<N;j++) {
                    tmp = *(a+k*N+j);
                    *(a+k*N+j) = *(a+u*N+j);
                    *(a+u*N+j)=tmp;
                }
            }

            // Now do block reduction
            *(p+k) = iPivot;
            if ( *(a+k*N+k) != ZERO ) {
		if(N < numThreads){
                	for (rows=k+1;rows<N;rows++) { 
                    		*(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);

                    		for (rows2=k+1;rows2<N;rows2++) { 
                        		*(a+rows*N+rows2) = *(a+rows*N+rows2) - *(a+rows*N+k) * *(a+k*N+rows2) ;
                    		}
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
                		thread_args = ( struct args * )  malloc(sizeof( struct args));
            			thread_args->N   = N;
				thread_args->k = k;
                		thread_args->startRow = startRow;
                		thread_args->stopRow = stopRow; 
               			thread_args->Aptr = a;

          		 	pthread_create( thread_id+i, NULL, &dls_thread_worker, thread_args );
            		    }
        		}
        		for(int i=0; i < numThreads ; i++) {
      	      			pthread_join( *(thread_id+i), NULL); 
        		}

        		free(numberOfRows);
        		free(thread_id);
		}
            }

            else {

                /* Handle the case of a zero pivot element, singular matrix */

                printf( "Element a[%d][%d} = %f\n", k, k, *(a+k*N+k)); 
                printf( " *** MATRIX A IS SINGULAR *** \n");
                printf( "    -- EXECUTION HALTED --\n");
                exit(1);
            }

        }
        // Now that we know we have reduced the matrices, start the 
        // back substitution process to solve for vector x.


        /* We now need to arrange b so that it has undergone the same 
         * operations as the matrix a.  This will then form
         * the vector y for the solution of Ux=y where U is the 
         * upper-triangular matrix formed in the elimination process
         * above. 
         */

	///if(N < numThreads){
        	for (k=0; k<N-1; k++ ) {
            		tmp = *(b+k);
            		*(b+k) = *(b+ *(p+k));
            		*(b+ *(p+k)) = tmp;

            		for (j=k+1;j<N;j++) 
                		*(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
        	} 
	/*}
	else{
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
                numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

                for (int i=0; i<numThreads; i++ ){
                	*(numberOfRows+i) = (N-1) / numThreads;
                }
                for (int i=0; i< (N-1) % numThreads; i++ ){
                        *(numberOfRows+i) = *(numberOfRows+i) + 1;
                }

                        stopRow=0;
                        for(int i=0; i < numThreads ; i++) {
                            {
                                startRow=stopRow;
                                stopRow=startRow+*(numberOfRows+i);
                                thread_args3 = ( struct args3 * )  malloc(sizeof( struct args3));
                                thread_args3->N   = N;
                                thread_args3->startRow = startRow;
                                thread_args3->stopRow = stopRow;
                                thread_args3->Aptr = a;
                                thread_args3->Bptr = b;
				thread_args3->Pptr = p;

                                pthread_create( thread_id+i, NULL, &dls_thread_worker3, thread_args3 );
                             }
                         }
                         for(int i=0; i < numThreads ; i++) {
                                pthread_join( *(thread_id+i), NULL);
                         }

                         free(numberOfRows);
                         free(thread_id);
	}*/
        // Now do the backward substitution to get the solution
        // vector x

        *(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
	//if(N < numThreads){
        	for (i=N-2;i>=0;i--){
            		tmp = 0.0;
            		for (j=i+1;j<N;j++) {
                		tmp = tmp + *(a+i*N+j) * *(b+j);
            		}
            		*(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
        	}
	/*}
	else{
                	thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
                        numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

                        for (int i=0; i<numThreads; i++ ){
                	        *(numberOfRows+i) = (N-2) / numThreads;
                        }
                        for (int i=0; i< (N-2) % numThreads; i++ ){
                                *(numberOfRows+i) = *(numberOfRows+i) + 1;
                        }

                        stopRow=0;
                        for(int i=0; i < numThreads ; i++) {
                            {
                                startRow=stopRow;
                                stopRow=startRow+*(numberOfRows+i);
                                thread_args2 = ( struct args2 * )  malloc(sizeof( struct args2));
                                thread_args2->N   = N;
                                thread_args2->startRow = startRow;
                                thread_args2->stopRow = stopRow;
                                thread_args2->Aptr = a;
				thread_args2->Bptr = b;

                         	pthread_create( thread_id+i, NULL, &dls_thread_worker2, thread_args2 );
                             }
                         }
                         for(int i=0; i < numThreads ; i++) {
                         	pthread_join( *(thread_id+i), NULL);
                         }

                         free(numberOfRows);
                         free(thread_id);
                   }*/

        //for (i=0;i<N;i++) *(x+i) = *(b+i);
        if(N < numThreads){
                for (i=0;i<N;i++) *(x+i) = *(b+i);
        }
        else{
                thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
                                        numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

                                        for (int i=0; i<numThreads; i++ ){
                                                *(numberOfRows+i) = (N) / numThreads;
                                        }
                                        for (int i=0; i< (N) % numThreads; i++ ){
                                                *(numberOfRows+i) = *(numberOfRows+i) + 1;
                                        }

                                        stopRow=0;
                                        for(int i=0; i < numThreads ; i++) {
                                           {
						//printf("startRow: %d\n", startRow);
                                                startRow=stopRow;
                                                stopRow=startRow+*(numberOfRows+i);
                                                thread_args5 = ( struct args5 *) malloc(sizeof( struct args5));
                                                thread_args5->startRow = startRow;
                                                thread_args5->stopRow = stopRow;
                                                thread_args5->Xptr = x;
                                                thread_args5->Bptr = b;

                                                pthread_create( thread_id+i, NULL, &dls_thread_worker6, thread_args5 );
                                           }
                                        }
                                        for(int i=0; i < numThreads ; i++) {
                                                pthread_join( *(thread_id+i), NULL);
                                        }

                                        free(numberOfRows);
                                        free(thread_id);
        }

        // At this point the solution to the system should be in vector x 

        free(p);
    }

    else {

        // Since we know the matrix is strictly diagonally dominant, verify
        // that none of the pivot elements are equal to zero

        singular = 1; 
        i=0;
        while ( i<N  && singular ) {
            singular = *(a+i*N+i) == ZERO;   
            i++;
        }

        if ( singular ) {
            printf( " *** MATRIX A IS SINGULAR *** \n");
            printf( "    -- EXECUTION HALTED -- \n");
            exit(1);
        }

        // We know at this point that we have a strictly diagonally dominant matrix that is
        // not singular -- so it sould be possible to do the LU factorization.
        //   (modified from Golub and van Loan, Chapter 3)


	//if(N < numThreads){
	        for (k=0; k<N-1; k++) {

		    if(N < numThreads){
        	    for (rows=k+1;rows<N;rows++) {
                	*(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);

      	            	for (rows2=k+1;rows2<N;rows2++) { 
                        	*(a+rows*N+rows2) = *(a+rows*N+rows2) - 
                        	*(a+rows*N+k) * *(a+k*N+rows2) ;
                    	}
                    }}
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
                                thread_args = ( struct args * )  malloc(sizeof( struct args));
                                thread_args->N   = N;
                                thread_args->k = k;
                                thread_args->startRow = startRow;
                                thread_args->stopRow = stopRow;
                                thread_args->Aptr = a;

                                pthread_create( thread_id+i, NULL, &dls_thread_worker, thread_args );
                            }
                        }
                        for(int i=0; i < numThreads ; i++) {
                                pthread_join( *(thread_id+i), NULL);
                        }

                        free(numberOfRows);
                        free(thread_id);
                }




                }
	/*}
	else{
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
                numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

                for (int i=0; i<numThreads; i++ ){
                	*(numberOfRows+i) = (N-1) / numThreads;
                }
                for (int i=0; i< (N-1) % numThreads; i++ ){
                                *(numberOfRows+i) = *(numberOfRows+i) + 1;
                }

                stopRow=0;
                        for(int i=0; i < numThreads ; i++) {
                           {
                                startRow=stopRow;
                                stopRow=startRow+*(numberOfRows+i);
                                thread_args = ( struct args * )  malloc(sizeof( struct args));
                                thread_args->N   = N;
                                thread_args->k = k;
                                thread_args->startRow = startRow;
                                thread_args->stopRow = stopRow;
                                thread_args->Aptr = a;

                                pthread_create( thread_id+i, NULL, &dls_thread_worker4, thread_args );
                            }
                        }
                        for(int i=0; i < numThreads ; i++) {
                                pthread_join( *(thread_id+i), NULL);
                        }

                        free(numberOfRows);
                 free(thread_id);

	}*/

        // At this point the LU factorizaton should be done and we have to do two
        // triangular back substitutions.  The solution to Ax=b is solved by first 
        // solving Ly=b for y and then Ux=y for the solution vector x.

        // Solving lower-triangular (Ly=b) first, overwriting b with y

        for (k=0; k<N-1; k++ ) {
            for (j=k+1;j<N;j++) 
                *(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
        } 

        // Now we can do the backward substitution to get the solution
        // vector x for the upper-triangular system (Ux=y) overwriting y (stored in b)
        // with x

        *(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
        for (i=N-2;i>=0;i--){
            tmp = 0.0;
            for (j=i+1;j<N;j++) {
                tmp = tmp + *(a+i*N+j) * *(b+j);
            }
            *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
        }


	if(N < numThreads){
        	for (i=0;i<N;i++) *(x+i) = *(b+i);
	}
	else{
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
					numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

					for (int i=0; i<numThreads; i++ ){
            					*(numberOfRows+i) = (N) / numThreads;
        				}
        				for (int i=0; i< (N) % numThreads; i++ ){
           					*(numberOfRows+i) = *(numberOfRows+i) + 1;
        				}

					stopRow=0;
        				for(int i=0; i < numThreads ; i++) {
            				   {   
						//printf("startRow: %d\n", startRow);
                				startRow=stopRow;
						//printf("startRow: %d\n", startRow);
                				stopRow=startRow+*(numberOfRows+i);
                				thread_args5 = ( struct args5 * )  malloc(sizeof( struct args5));
                				thread_args5->startRow = startRow;
                				thread_args5->stopRow = stopRow; 
               			 		thread_args5->Xptr = x;
						thread_args5->Bptr = b;

            			    		pthread_create( thread_id+i, NULL, &dls_thread_worker6, thread_args5 );
            		    		   }
        				}
        				for(int i=0; i < numThreads ; i++) {
      			      			pthread_join( *(thread_id+i), NULL); 
        				}

        				free(numberOfRows);
        				free(thread_id);
	}
        // At this point the solution to the system should be in vector x 


    }

}


int strictlyDiagonallyDominant( int N, double *a ) {

    double sum;
    int i, testPassed, row;

    testPassed = 1;
    row = 0;
    sum = 0.0;
    for (row=0;row<N;row++) { 
        if (testPassed) {
            sum = 0.0;
            for (i=0;i<N;i++) sum+=*(a+row*N+i);
            sum-=fabs(*(a+row*N+row)); 
            testPassed = fabs(*(a+row*N+row)) > sum;
        }
    }

    return testPassed;
}


void *dls_thread_worker( struct args *thread_args  ) {

    int i, j, k;
    double val, tmp, sum;
    int rows, rows2, testPassed, rowStart, rowStop, N;
    double *a;

    N        =  thread_args->N;
    k = thread_args->k;
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow;
    a        =  thread_args->Aptr;

    //printf("k: %d rowStart: %d rowStop: %d\n", k, rowStart, rowStop);

    for (rows=rowStart;rows<rowStop;rows++) {
    	*(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);
        for (rows2=k+1;rows2<N;rows2++) {
        	*(a+rows*N+rows2) = *(a+rows*N+rows2) - *(a+rows*N+k) * *(a+k*N+rows2) ;
        }
    }


    free(thread_args);
    pthread_exit(NULL);
}

void *dls_thread_worker2( struct args2 *thread_args  ) {

    int i, j, k;
    double val, tmp, sum;
    int rows, rows2, testPassed, rowStart, rowStop, N;
    double *a, *b;

    N        =  thread_args->N;
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow;
    a        =  thread_args->Aptr;
    b = thread_args->Bptr;

    //printf("k: %d rowStart: %d rowStop: %d\n", k, rowStart, rowStop);

    for (i=rowStop;i>=rowStart;i--){
            tmp = 0.0;
            for (j=i+1;j<N;j++) {
                tmp = tmp + *(a+i*N+j) * *(b+j);
            }
            *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
     }

    free(thread_args);
    pthread_exit(NULL);
}

void *dls_thread_worker3( struct args3 *thread_args  ) {

    int i, j, k;
    double val, tmp, sum;
    int rows, rows2, testPassed, rowStart, rowStop, N;
    double *a, *b;
    int *p;

    N        =  thread_args->N;
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow;
    a        =  thread_args->Aptr;
    b = thread_args->Bptr;
    p = thread_args->Pptr;

    //printf("k: %d rowStart: %d rowStop: %d\n", k, rowStart, rowStop);

    for (k=rowStart; k<rowStop; k++ ) {
    	tmp = *(b+k);
        *(b+k) = *(b+ *(p+k));
        *(b+ *(p+k)) = tmp;

        for (j=k+1;j<N;j++)
        	*(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);
    }


    free(thread_args);
    pthread_exit(NULL);
}

void *dls_thread_worker4( struct args *thread_args  ) {

    int i, j, k;
    double val, tmp, sum;
    int rows, rows2, testPassed, rowStart, rowStop, N;
    double *a, *b;

    N        =  thread_args->N;
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow;
    a        =  thread_args->Aptr;

    //printf("k: %d rowStart: %d rowStop: %d\n", k, rowStart, rowStop);

	for (k=0; k<N-1; k++) {
                    for (rows=k+1;rows<N;rows++) {
                        *(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);

                        for (rows2=k+1;rows2<N;rows2++) {
                                *(a+rows*N+rows2) = *(a+rows*N+rows2) -
                                *(a+rows*N+k) * *(a+k*N+rows2) ;
                        }
                    }
        }
    free(thread_args);
    pthread_exit(NULL);

}

void *dls_thread_worker5( struct args4 *thread_args  ) {

    int i, j, k;
    double val, tmp, sum;
    int rows, rows2, testPassed, rowStart, rowStop, N;
    double *a, *b;
    int *p;

    N        =  thread_args->N;
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow;
    p = thread_args->Pptr;

    //printf("k: %d rowStart: %d rowStop: %d\n", k, rowStart, rowStop);

    for(k=rowStart;k<rowStop;k++) *(p+k)=k;

    free(thread_args);
    pthread_exit(NULL);
}

void *dls_thread_worker6( struct args5 *thread_args  ) {

    int i, rowStop, rowStart;
    double *x, *b;

    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow;
    x = thread_args->Xptr;
    b = thread_args->Bptr;

    //printf("rowStart: %d rowStop: %d\n", rowStart, rowStop);

    for (i=rowStart;i<rowStop;i++) *(x+i) = *(b+i);    

    free(thread_args);
    pthread_exit(NULL);
}


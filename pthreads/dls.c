#ifdef __cplusplus
extern "C" {
#endif
	void dls_(int *num_threads, int *len, double* mat, double* vec, double* rvec);
#ifdef __cplusplus
}
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

void *dls_thread_worker();

struct args {
    int N; 
    int startRow;
    int stopRow;
    double *aptr;
    double *bptr;
    double *cptr;
};

int strictlyDiagonallyDominant(int N, double *a);

void dls_(int* num_threads, int* len, double* a, double* b, double* x){

	int numThreads = *num_threads;
	int matrixDimension = *len;
	int *numberOfRows;
	int startRow, stopRow;
	pthread_t *thread_id;
	struct args *thread_args;

	if(matrixDimension < numThreads){
		int i, j, k, N, u;
		int singular, iPivot, rows, rows2;
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
					if(fabs(*(a+u*N+k)) > fabs(pivotMax)){
						pivotMax = *(a+u*N+k);
						iPivot = u;
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
							*(a+rows*N+rows2) = *(a+rows*N+rows2) - *(a+rows*N+k) * *(a+k*N+rows2) ;
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
				*(b+ *(p+k)) = tmp;
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

	else{

		thread_id = (pthread_t *) malloc(numThreads * sizeof(pthread_t));
		numberOfRows = (int*)malloc(numThreads*sizeof(int));

		for(int i = 0; i<numThreads;i++){
			*(numberOfRows+i) = matrixDimension/numThreads;
		}

		for(int i =0; i<matrixDimension%numThreads;i++){
			*(numberOfRows+i) = *(numberOfRows+i) + 1;
		}
		
		stopRow = 0;

		for(int i = 0; i<numThreads; i++){
			{
			startRow = stopRow;
			stopRow=startRow+*(numberOfRows+i);
			thread_args = ( struct args * )  malloc(sizeof( struct args));
                	thread_args->N   = matrixDimension;
                	thread_args->startRow = startRow;
                	thread_args->stopRow = stopRow; 
                	thread_args->aptr = a;
                	thread_args->bptr = b;
                	thread_args->cptr = c;

                	pthread_create( thread_id+i, NULL, &dls_thread_worker, thread_args );
            		}
		}
		for(int i = 0; i <numThreads; i++){
			pthread_join( *(thread_id+i), NULL);
		}

		free(numberOfRows);
		free(thread_id);
	}

}	

void *dls_thread_worker(struct args *thread_args){
	int i, j, k;
    	double val;
    	int rowStart, rowStop, N; 
    	double *a, *b, *c;

	N = thread_args->N;
    	rowStart = thread_args->startRow;
    	rowStop = thread_args->stopRow; 
    	a = thread_args->aptr;
    	b = thread_args->bptr;
    	c = thread_args->cptr;

	int u;
	int singular, iPivot, rows, rows2;
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
				if(fabs(*(a+u*N+k)) > fabs(pivotMax)){
					pivotMax = *(a+u*N+k);
					iPivot = u;
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
						*(a+rows*N+rows2) = *(a+rows*N+rows2) - *(a+rows*N+k) * *(a+k*N+rows2) ;
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
			*(b+ *(p+k)) = tmp;
			for (j=k+1; j<N; j++){
				*(b+j) = *(b+j) - *(b+k) * *(a+N*j+k);
			}
		}

		*(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
		for (i=rowStop-2;i>=rowStart;i--){
			tmp = 0.0;
			for(j=i+1;j<N;j++){
				tmp = tmp + *(a+i*N+j) * *(b+j);
			}
			*(b+i) = (*(b+i) - tmp)/ *(a+i*N+i);
		}

		for (i=rowStart;i<rowStop;i++){
			*(x+i) = *(b+i);
		}

		free(p);
	}

	else {
		singular = 1;
		i=rowStart;
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
        	for (i=rowStop-2;i>=rowStart;i--){
            		tmp = 0.0;
            		for (j=i+1;j<N;j++) {
                		tmp = tmp + *(a+i*N+j) * *(b+j);
            		}
            		*(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
        	}

        	for (i=rowStart;i<rowStop;i++) *(x+i) = *(b+i);

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


#ifdef __cplusplus
extern "C" {
#endif
	void dot_(int *threads, int *N, double *vec1, double *vec2, double *rvec);
#ifdef __cplusplus
	}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

void *dot_thread_worker();

struct args {
	int startRow;
	int stopRow;
	double *v1ptr;
	double *v2ptr;
	double *rvalptr;
};

void dot_(int *threads, int *N, double *vec1, double *vec2, double *rval){

	int numThreads = *threads;
	int len = *N;
	int *numberOfRows;
	int startRow, stopRow;
	pthread_t *thread_id;
	struct args *thread_args;


	int i;
	if (len < numThreads){
		
        	for(i=0; i < len; i++){

                	*rval += (*(vec1+i) * *(vec2+i));
        
		}
	}	
	else{
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
		numberOfRows = (int *) malloc (numThreads * sizeof(int));
		for(int i = 0; i < numThreads; i++){
			*(numberOfRows+i) = len/numThreads;
		}
		
		for(int i = 0; i < len % numThreads; i++){
			*(numberOfRows+i) = *(numberOfRows+i) + 1;
		}

		stopRow = 0;
		for (int i = 0; i < numThreads; i++){
			startRow = stopRow;
			stopRow = startRow + *(numberOfRows+i);
			thread_args = (struct args *) malloc(sizeof(struct args));
			thread_args->startRow = startRow;
			thread_args->stopRow = stopRow;
			thread_args->v1ptr = vec1;
			thread_args->v2ptr = vec2;
			thread_args->rvalptr = rval;

			pthread_create(thread_id+i, NULL, &dot_thread_worker, thread_args);
		}
	}
	for(int i = 0; i < numThreads; i++){
		pthread_join( *(thread_id+i), NULL);
	}

	free(numberOfRows);
	free(thread_id);
}

void *dot_thread_worker(struct args *thread_args){

	int i, j, k;
	double val;
	int rowStart, rowStop, N;
	double *vec1, *vec2, *rval;

	rowStart = thread_args->startRow;
	rowStop = thread_args->stopRow;
	vec1 = thread_args->v1ptr;
	vec2 = thread_args->v2ptr;
	rval = thread_args->rvalptr;

	for(i=rowStart; i < rowStop; i++){
        	*rval += (*(vec1+i) * *(vec2+i));
        }
	
	free(thread_args);
	pthread_exit(NULL);
}

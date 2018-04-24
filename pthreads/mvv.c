#ifdef __cplusplus
extern "C" {
#endif
	void mvv_(int *threads, int *N, double* mat, double* vec, double* vresults);
#ifdef __cplusplus
	}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

struct args {
	int N;
	int startRow;
	int stopRow;
	double *matptr;
	double *vecptr;
	double *vresultsptr;
};


void *mvv_thread_worker(struct args *thread_args);

void mvv_(int *threads, int *N, double* mat, double* vec, double* vresults){

	int i, j;

	int numThreads = *threads;

        //printf("numthreads: %d\n", numThreads);

	//printf("", );

	int length = *N;

	//printf("length: %d\n", length);

	int *numberOfRows;
	int startRow, stopRow;
	pthread_t *thread_id;
	struct args *thread_args;



	// single threaded
	if( length < numThreads){
		
		for(i=0;i<length;i++){
	                for(j=0;j<length;j++){
        	                *(vresults+i) += *(mat+(length*i)+j) * *(vec+j);
                	        //printf("%f\n", *(vresults+i));
                	}
        	}
	}

	else{
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
                numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

                for (int i=0; i<numThreads; i++ ){
                              *(numberOfRows+i) = length / numThreads;
                }
                for (int i=0; i< length % numThreads; i++ ){
                              *(numberOfRows+i) = *(numberOfRows+i) + 1;
                }

                stopRow=0;
                for(int i=0; i < numThreads ; i++) {
                    {
                	startRow=stopRow;
                        stopRow=startRow+*(numberOfRows+i);
                        thread_args = ( struct args * )  malloc(sizeof( struct args));
			thread_args->N = length;
                        thread_args->startRow = startRow;
                        thread_args->stopRow = stopRow;
                        thread_args->matptr = mat; 
			thread_args->vecptr = vec;
		 	thread_args->vresultsptr = vresults;

                        pthread_create( thread_id+i, NULL, &mvv_thread_worker, thread_args );
                     }
                 }
                 for(int i=0; i < numThreads ; i++) {
                 	pthread_join( *(thread_id+i), NULL);
                 }

                 free(numberOfRows);
                 free(thread_id);		
	}


}

void *mvv_thread_worker(struct args *thread_args){
        int i, j, k;
        double val;
        int rowStart, rowStop, N;
        double *mat, *vec, *vresults;

        N = thread_args->N;
        rowStart = thread_args->startRow;
        rowStop = thread_args->stopRow;
        mat = thread_args->matptr;
        vec = thread_args->vecptr;
        vresults = thread_args->vresultsptr;

	//printf("N: %d rowStart: %d rowStop: %d\n", N, rowStart, rowStop);
	/*
        for(i=rowStart;i<rowStop;i++){
                for(j=0;j<N;j++){
                        *(vresults+i) += *(mat+(N*i)+j) * *(vresults+j);
                }
        }*/

	for(i=rowStart;i<rowStop;i++){
        	for(j=0;j<N;j++){
                	*(vresults+i) += *(mat+(N*i)+j) * *(vec+j);
                        //printf("%f\n", *(vresults+i));
                }
        }



	free(thread_args);
	pthread_exit(NULL);
}



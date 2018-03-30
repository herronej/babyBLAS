#ifdef __cplusplus
extern "C" {
#endif
	void mvv_(int *thread, int *N, double* mat, double* vec, double* vresults);
#ifdef __cplusplus
	}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

void *mvv_thread_worker();

struct args {
	int N;
	int startRow;
	int stopRow;
	double *matptr;
	double *vecptr;
	double *vresultsptr;
};

void mvv_(int *N, double* mat, double* vec, double* vresults){

	int numThreads = *threads;
	int matrixDimension = *len;
	int *numberOfRows;
	int startRow, stopRow;
	pthread_t *thread_id;
	struct args *thread_args;

	// single threaded
	if( matrixDimension < numThreads){
		for(int i = 0; i < matrixDimension; i++){
			for(int j = 0; j < matrixDimension; j++){
				*(vresults+i) += *(mat+(length*i)+j) * *(vresults+j);
			}
		}
	}

	else{
		thread_id = (pthread_t *) malloc(numThreads * sizeof(pthread_t));
		numberOfRows = (int*) malloc(numThreads * sizeof(int));

		for(int i = 0; i < numThreads; i++){
			*(numberOfRows+i) = matrixDimension / numThreads;
		}
		for(int i = 0; i < matrixDimension % numThreads; i++){
			*(numberOfRows+i) = *(numberOfRows+i) + 1;
		}

		stopRow = 0;

		for(int i = 0; i < numThreads; i++){
			{
			startRow = stopRow;
			stopRow = startRow+*(numberOfRows+i);
			thread_args = (struct args *) malloc(sizeof(struct args));
			thread_args->N = matrixDimension;
			thread_args->startRow = startRow;
			thread-args->stopRow = stopRow;
			thread_args->matptr = mat;
			thread_args->vecptr = vec;
			thread_args->vresultsptr = vresults;

			pthread_create( thread_id+i, NULL, &mvv_thread_worker, thread_args);}
		}

		for(int i = 0; i < numThreads; i++){
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

        for(i=rowStart;i<rowStop;i++){
                for(j=0;j<N;j++){
                        *(vresults+i) += *(mat+(length*i)+j) * *(vresults+j);
                }
        }

	free(thread_args);
	pthread_exit(NULL);
}



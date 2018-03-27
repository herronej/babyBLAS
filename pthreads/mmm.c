#ifdef __cplusplus
extern "C" {
#endif
	void mmm_(int *threads, int *len, double *a, double *b, double *c);
#ifdef __cplusplus
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

void *mmm_thread_worker();

struct args {
	int N;
	int startRow;
	int stopRow;
	double *Aptr;
	double *Bptr;
	double *Cptr;
};

void mmm_(int *threads, int *len, double *A, double *B, double *C){
	
	int numThreads = *threads;
	int matrixDimension = *len;
	int *numberOfRows;
	int startRow, stopRow;
	pthread_t *thread_id;
	struct args *thread_args;

	// do single matrix multiplication if fewer dimensions than threads
	if(matrixDimension < numThreads){
		for(int i = 0; i < matrixDimension; i++){
			for(int j = 0; j < matrixDimension; j++){
				*(C+(i*matrixDimension+j))=0.0;
				for (int k = 0; k < matrixDimension;k++){
					*(C+(i*matrixDimension+j)) += *(A+(i*matrixDimension+k)) * *(B+(k*matrixDimension+j));
				}
			}
		}
	}
	else{

		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));

		numberOfRows = (int *) malloc(numThreads*sizeof(int));
		
		for(int i = 0; i < numThreads; i++){
			*(numberOfRows+i) = matrixDimension /numThreads;
		}

		for(int i = 0; i < matrixDimension % numThreads; i++){
			*(numberOfRows+i) = *(numberOfRows+i) + 1;
		}


		stopRow = 0;
		for(int i = 0; i < numThreads; i++){
			startRow = stopRow;
			stopRow = startRow + *(numberOfRows+i);
			thread_args = (struct args *) malloc(sizeof(struct args));
			thread_args->N = matrixDimension;
			thread_args->startRow = startRow;
			thread_args->stopRow = stopRow;
			thread_args->Aptr = A;
			thread_args->Bptr = B;
			thread_args->Cptr = C;

			pthread_create(thread_id+i, NULL, &mmm_thread_worker, thread_args);
		

		}
	}
	for(int i = 0; numThreads; i++){
		pthread_join(*(thread_id+i), NULL);
	}
	free(numberOfRows);
	free(thread_id);

    }
}

void *mmm_thread_worker(struct args *thread_args){
	
	int i, j, k;
	double val;
	int rowStart, rowStop, N;
	double *A, *B, *C;

	N = thread_args->N;
	rowStart = thread_args->startRow;
	rowStop = thread_args->stopRow;
	A = thread_args->Aptr;
	B = thread_args->Bptr;
	C = thread_args->Cptr;

	for(i=rowStart;i<rowStop;i++){
		for(j=0; j<N; j++){
			*(C+(i*N+j))=0.0;
			for(k=0; k < N; k++){
				*(C+(i*N+j)) += *(A+(i*N+k)) * *(B+(k*N+j));
			}
		}
	}

	free(thread_args);
	pthread_exit(NULL);

}
	

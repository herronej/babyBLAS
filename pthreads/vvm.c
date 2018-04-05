#ifdef __cplusplus
extern "C"{
#endif
	void vvm_(int *threads, int *len, double *va, double *vb, double *ma);
#ifdef __cplusplus
	}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


void *vvm_thread_worker();

struct args {
	int N;
	int startRow;
	int stopRow;
	double *vaptr;
	double *vbptr;
	double *maptr;	
};

void vvm_(int *threads, int *len, double *va, double *vb, double *ma){

	int numThreads = *threads;
	int matrixDimension = *len;
	int *numberOfRows;
	int startRow, stopRow;
	pthread_t *thread_id;
	struct args *thread_args;

	
	int alength = *len;

	//int i, j;

	if(matrixDimension < numThreads){
		for(int i=0;i<alength;i++){
			for(int j=0;j<alength;j++){
				*(ma+(alength*i)+j) = *(va+i) * *(vb+j);
			}
		}
	}
	else{
		thread_id = (pthread_t *) malloc(numThreads * sizeof(pthread_t));
		numberOfRows = (int*) malloc(numThreads*sizeof(int));
		
		for(int i = 0; i < numThreads; i++){
			*(numberOfRows+i) = matrixDimension / numThreads;
		}

		for(int i = 0; i < matrixDimension % numThreads; i++){
			*(numberOfRows+i) = *(numberOfRows+i) + 1;
		}

		stopRow = 0;

		for(int i = 0; i < numThreads; i++){
			{startRow=stopRow;
			stopRow=startRow+*(numberOfRows+i);
			thread_args = (struct args*) malloc(sizeof(struct args));
			thread_args->N = matrixDimension;
			thread_args->startRow = startRow;
			thread_args->stopRow = stopRow;
			thread_args->vaptr = va;
			thread_args->vbptr = vb;
			thread_args->maptr = ma;
			pthread_create(thread_id+i, NULL, &vvm_thread_worker, thread_args);
			}
		}
		for(int i = 0; i < numThreads; i++){
			pthread_join( *(thread_id+i), NULL);
		}

		free(numberOfRows);
		free(thread_id);
	}

}

void *vvm_thread_worker(struct args *thread_args){
	
	//int i, j;
	double val;
	int rowStart, rowStop, N;
	double *va, *vb, *ma;

	N = thread_args->N;
	rowStart = thread_args->startRow;
	rowStop = thread_args->stopRow;
	va = thread_args->vaptr;
	vb = thread_args->vbptr;
	ma = thread_args->maptr; 

	for(int i=rowStart;i<rowStop;i++){
        	for(int j=0;j<N;j++){
                        *(ma+(N*i)+j) = *(va+i) * *(vb+j);
                }
        }

	free(thread_args);
	pthread_exit(NULL);
}

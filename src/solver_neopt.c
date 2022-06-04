/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	int i, j, k;
	
	double *C = (double *) calloc(N*N, sizeof(double));
	double *aux1 = (double *) calloc(N*N, sizeof(double));
	double *aux2 = (double *) calloc(N*N, sizeof(double));
	double *aux3 = (double *) calloc(N*N, sizeof(double));
	double *At = (double *) calloc(N*N, sizeof(double));
	double *Bt = (double *) calloc(N*N, sizeof(double));

	// Calculate A transpose
	for (i = 0; i < N; i++) {
		for (j = i; j < N; j++) {
			At[j * N + i] = A[i * N + j];
		}
	}

	// Calculate B transpose
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			Bt[j * N + i] = B[i * N + j];
		}
	}

	// aux1 = A * At
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = i; k < N; k++) {
				aux1[i * N + j] += A[i * N + k] * At[k * N + j];	
			}
		}
	}

	// aux2 = Bt * B
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				aux2[i * N + j] += Bt[i * N + k] * B[k * N + j];
			}
		}
	}

	// aux3 = B * aux1
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				aux3[i * N + j] += B[i * N + k] * aux1[k * N + j]; 
			}
		}
	}	

	// C = aux3 + aux2
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			C[i * N + j] += aux3[i * N + j] + aux2[i * N + j];
		}
	}

	free(aux1);
	free(aux2);
	free(aux3);
	free(At);
	free(Bt);

	return C;
}

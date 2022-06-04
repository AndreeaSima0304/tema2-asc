/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	register int i = 0;
       	register int j = 0;
	register int k = 0;

	double *C = (double *) calloc(N*N, sizeof(double));
        register double *aux1 = (double *) calloc(N*N, sizeof(double));
        register double *aux2 = (double *) calloc(N*N, sizeof(double));
        register double *aux3 = (double *) calloc(N*N, sizeof(double));
	register double *At = (double *) calloc(N*N, sizeof(double));
        register double *Bt = (double *) calloc(N*N, sizeof(double));
	
	// Calculate A transpose
	for (i = 0; i < N; ++i) {
		for (j = i; j < N; ++j) {
			At[j * N + i] = A[i * N + j];
		}
	}

	// Calculate B transpose
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			Bt[j * N + i] = B[i * N + j];
		}
	}

	// aux1 = A * At
	for (i = 0; i < N; ++i) { 
		register double *orig_pa = &A[i * N + i];
		for (j = 0; j < N; ++j) {
			register double *pa = orig_pa;
			register double *pat = &At[i * N + j];
			register double sum = 0;
			for (k = i; k < N; ++k) { 	
				sum += *pa * *pat;
				pa++;
				pat += N;
			}
			aux1[i * N + j] = sum;
                }
	}
	free(At);

	// aux2 = Bt * B
	for (i = 0; i < N; ++i) {
		register double *orig_pbt = &Bt[i * N + 0];
		for (j = 0; j < N; ++j) {
			register double *pbt = orig_pbt;
			register double *pb = &B[0 * N + j];
			register double sum = 0;
			for (k = 0; k < N; ++k) {
				sum += *pbt * *pb;
				pbt++;
				pb += N;
			}
			aux2[i * N + j] = sum;
		}
	}
	free(Bt);

	// aux3 = B * aux1
	for (i = 0; i < N; ++i) {
		register double *orig_pb = &B[i * N + 0];
		for (j = 0; j < N; ++j) {
			register double *pb = orig_pb;
			register double *paux1 = &aux1[0 * N + j];
			register double sum = 0;
			for (k = 0; k < N; ++k) {
				sum += *pb * *paux1;
				pb++;
				paux1 += N;
			}
			aux3[i * N + j] = sum;
		}
	}
	free(aux1);

	// C = aux3 + aux2
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			C[i * N + j] += aux3[i * N + j] + aux2[i * N + j];
		}
	}

	free(aux2);
	free(aux3);

	return C;
}

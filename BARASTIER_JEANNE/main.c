#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>
#include "param.h"

#define PHI_MASK ((uint64_t)~(uint64_t)0) // 2^64-1

int64_t* genpoly(int64_t n, int64_t rho) {
	int64_t i;
	int j;
	int64_t* poly = (int64_t*)malloc(n*sizeof(int64_t));

	/* J'ai vu que la fonction rand() ne génère en général que des nombres sur 15 bits,
	donc que vais faire appel à rand() 4 fois et concaténer les résultats afin d'avoir
	un nombre sur 60 bits (car rho peut être un très grand nombre) */
	for(i=0;i<n;i++){
		int64_t coeff;
		do {
			coeff = 0;
			for(j=0;j<4;j++){
				coeff = (coeff << 15)^(rand() & 0x7FFF); // génère un nombre, garde seulement les 15 bits de poids faible et concatène
			}
			coeff = coeff % (2*rho + 1) - rho; // replace mon coeff entre -rho et rho (exclusivement)
		} while (coeff <= -rho || coeff >= rho);
		poly[i] = coeff;
	}
	return poly;
}

__int128* poly_mul_modE_int128(const int64_t* A, const int64_t* B, const int64_t* E, int64_t n) {
	size_t tlen = 2*n - 1;
	__int128* T = (__int128*)calloc(tlen, sizeof(__int128));

	for (int64_t i = 0; i < n; i++) {
		for (int64_t j = 0; j < n; j++) {
			T[i+j] += (__int128)A[i] * (__int128)B[j];
		}
	}

	// réduction modulo E(X)
	for (int64_t i = 2*n - 2; i >= n; i--) {
		if (T[i] != 0) {
			__int128 coef = T[i];
			for (int64_t j = 0; j < n; j++) {
				T[i - n + j] -= coef * (__int128)E[j];
			}
			T[i] = 0;
		}
	}

	__int128* C = (__int128*)malloc(n * sizeof(__int128));
	for (int64_t i = 0; i < n; i++) C[i] = T[i];

	free(T);
	return C;
}

uint64_t* poly_mul_modE_then_modphi(const __int128* C128, const uint64_t* Minv, const int64_t* E, int64_t n) {
	size_t tlen = 2*n - 1;
	__int128* T = (__int128*)calloc(tlen, sizeof(__int128));

	for (int64_t i = 0; i < n; i++) {
		for (int64_t j = 0; j < n; j++) {
			unsigned __int128 uC = (unsigned __int128)C128[i];
			unsigned __int128 uM = (unsigned __int128)Minv[j];
			unsigned __int128 prod = uC * uM;
			T[i+j] += (__int128)prod;
		}
	}

	// reduction modulo E
	for (int64_t i = 2*n - 2; i >= n; i--) {
		if (T[i] != 0) {
			__int128 coef = T[i];
			for (int64_t j = 0; j < n; j++) {
				T[i - n + j] -= coef * (__int128)E[j];
			}
			T[i] = 0;
		}
	}

	uint64_t* Q64 = (uint64_t*)malloc(n * sizeof(uint64_t));
	for (int64_t i = 0; i < n; i++) {
		unsigned __int128 u = (unsigned __int128)T[i];
		Q64[i] = (uint64_t)(u & PHI_MASK);
	}

	free(T);
	return Q64;
}

__int128* poly_mul_modE_from_QM(const uint64_t* Q64, const int64_t* M, const int64_t* E, int64_t n) {
	size_t tlen = 2*n - 1;
	__int128* T = (__int128*)calloc(tlen, sizeof(__int128));

	for (int64_t i = 0; i < n; i++) {
		for (int64_t j = 0; j < n; j++) {
			unsigned __int128 uQ = (unsigned __int128)Q64[i];
			__int128 prod = (__int128)uQ * (__int128)M[j];
			T[i+j] += prod;
		}
	}

	// reduction mod E
	for (int64_t i = 2*n - 2; i >= n; i--) {
		if (T[i] != 0) {
			__int128 coef = T[i];
			for (int64_t j = 0; j < n; j++) {
				T[i - n + j] -= coef * (__int128)E[j];
			}
			T[i] = 0;
		}
	}

	__int128* S128 = (__int128*)malloc(n * sizeof(__int128));
	for (int64_t i = 0; i < n; i++) S128[i] = T[i];

	free(T);
	return S128;
}

int64_t* MONTGMULT(const int64_t* A, const int64_t* B, const int64_t* E,
				   const int64_t* M, const uint64_t* Minv, int64_t n)
{
	// C128 = A*B mod E  (in __int128)
	__int128* C128 = poly_mul_modE_int128(A, B, E, n);

	// Q64 = (C * Minv mod E) mod 2^64
	uint64_t* Q64 = poly_mul_modE_then_modphi(C128, Minv, E, n);

	// S128 = Q * M mod E
	__int128* S128 = poly_mul_modE_from_QM(Q64, M, E, n);

	// R[i] = (C128[i] - S128[i]) / 2^64  (exact division in PMNS)
	int64_t* R = (int64_t*)malloc(n * sizeof(int64_t));

	for (int64_t i = 0; i < n; i++) {
		__int128 numer = C128[i] - S128[i];
		R[i] = (int64_t)(numer >> 64);
	}

	free(C128);
	free(Q64);
	free(S128);
	return R;
}

// pour m'aider à print
void print_poly_int64(const int64_t* P, int64_t n) {
	for (int64_t i = 0; i < n; i++) {
		printf("%" PRId64 "%s", P[i], (i < n-1) ? ", " : "");
	}
}

int main(void) {
	srand((unsigned)time(NULL));

	for (int iter = 0; iter < 100; iter++) {
		int64_t* A = genpoly(n, rho);
		int64_t* B = genpoly(n, rho);
		int64_t* Cprime = MONTGMULT(A, B, E, M, Minv, n);

		printf("\nItération %d\n", iter + 1);
		printf("A(X) = ");
		print_poly_int64(A, n);
		printf("\nB(X) = ");
		print_poly_int64(B, n);
		printf("\nC'(X) = ");
		print_poly_int64(Cprime, n);
		printf("\n");

		free(A);
		free(B);
		free(Cprime);
	}

	return 0;
}
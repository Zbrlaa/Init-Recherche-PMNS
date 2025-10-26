#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> //Pour utiliser int64_t
#include <time.h>
#include "param.h"


void genpoly(int64_t *poly, int n) {
	for(int i = 0; i < n; i++) {
		// Generate coefficients in [-rho+1, rho-1]
		poly[i] = (rand() % (2 * rho - 1)) - (rho - 1);
	}
}


void print_poly(const int64_t *poly, int n) {
	for(int i = 0; i < n; i++) {
		if(i == 0) printf("%ld", poly[i]);
		else printf(" + %ldX^%d", poly[i], i);
	}
	printf("\n");
}

void MONTGMUL(int64_t *Cprime, const int64_t *A, const int64_t *B) {
	//var temporaire, pas de malloc, plus rapide
	__int128_t C[2*n-1];
	__int128_t Q[n];
	__int128_t temp[n];

	//Init
	for(int i = 0; i < 2*n-1; i++) C[i] = 0;
	for(int i = 0; i < n; i++) Q[i] = temp[i] = 0;

	//C = A × B
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			C[i+j] += (__int128_t)A[i] * (__int128_t)B[j];
		}
	}

	//Q = (C × Minv) % PHI
	for(int i = 0; i < n; i++) {
		__int128_t sum = 0;
		for(int j = 0; j < n; j++) {
			sum += C[j] * (__int128_t)Minv[i];
		}
		Q[i] = (int64_t)(sum & PHI_MASK);
	}

	//temp = Q × M
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			temp[i] += Q[j] * (__int128_t)M[i];
		}
	}

	//C' = (C - temp) / PHI
	for(int i = 0; i < n; i++) {
		Cprime[i] = (int64_t)((C[i] - temp[i]) >> 64);
	}
}


int main() {
	srand(time(NULL));
	
	int64_t A[n], B[n], Cprime[n];
	
	for(int i = 0; i < 100; i++){
		printf("Itération n°%d :\n", i);

		genpoly(A, n);
		genpoly(B, n);
		
		printf("A = ");
		print_poly(A, n);
		printf("B = ");
		print_poly(B, n);
		
		MONTGMUL(Cprime, A, B);
		
		printf("C' = ");
		print_poly(Cprime, n);

		printf("\n");
	}
	
	return 0;
}
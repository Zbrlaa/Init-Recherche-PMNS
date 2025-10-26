from sage.all import *
from random import *
	

def sqandmult(a,b,p):
	c = 1
	while b :
		if b&1 :
			c = (c*a)%p
		a = (a*a)%p
		b = b >> 1

	return c


def pmns(PSIZE):
	p = random_prime(2**PSIZE)
	n = floor(log(p,2)/64)+1
	phi = 2**64
	lambd_max = 200

	while(True) :
		print(f"n : {n}")
		lambd = 2

		while(2*n*abs(lambd)*p**(1/n) < phi and lambd < lambd_max) :
			K = GF(p)
			pol = PolynomialRing(K, "X")
			X = pol("X")
			E = ZZ["X"](f"X^{n}-{lambd}")
			Xp = sqandmult(X,p,E)
			d = xgcd((Xp-X)%E, E)[0]

			#pas de fact commun
			if d == 1 or d[0]==0 :
				lambd += 2
				continue

			gamma = int(-d[0]%p)
			
			Bmat = matrix(ZZ,n,n)
			Bmat[0,0] = p
			for i in range(1,n):
				Bmat[i,0] = -gamma**i
				Bmat[i,i] = 1
			B = Bmat.LLL()

			Bn = B.norm(1)
			if 2*n*abs(lambd)*Bn < phi :
				rho = int(Bn - 1)
				print(f"Rho : {rho}")

			#Echec calcul rho
			else :
				lambd += 2
				continue

			M = None
			for T in B:
				if T[0] % 2:
					M = T
					break

			M = ZZ["X"](list(M))
			print(f"M(X)={M}")

			dM, uM, _ = xgcd(M, E)
			dM = int(dM)
			dm1 = pow(dM, -1, phi)

			Mm1 = dm1 * uM%phi
			Mm1 = ZZ["X"](Mm1)
			print(f"M-1(X)={Mm1}")

			coeffs_A = [randrange(-rho+1,rho) for _ in range(n)]
			coeffs_B = [randrange(-rho+1,rho) for _ in range(n)]
			A = ZZ["X"](coeffs_A)
			B = ZZ["X"](coeffs_B)
			# print(f"A(X)={A}")
			# print(f"B(X)={B}")

			gamma = int(gamma)

			C = (A*B)%E
			# print(f"C(X)={C}")

			Q = ((C*Mm1)%E)%phi
			# print(f"Q(X)={Q}")

			Cp = (C - (Q*M)%E)/phi
			Cp = ZZ["X"](Cp)
			# print(f"C'(X)={Cp}")

			print(f"E(X) = {E}")
			result = (Cp(gamma)%p == (A(gamma)*B(gamma)*pow(phi,-1,p))%p)
			print(f"Verif pour n({n}) et lambda({lambd}): {result}")

			return
		
		n += 1

pmns(PSIZE = 200)
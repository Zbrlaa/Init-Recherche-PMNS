from sage.all import *
import random

def sqandmult(a,b,p):
	c = 1
	while b:
		if b&1:
			c = (c*a)%p
		a = (a*a)%p
		b = b>>1
	return c

def matrice_base_courte(p, n, gamma):
	Bmat = matrix(ZZ,n,n)
	Bmat[0,0] = p
	for i in range(1,n):
		Bmat[i,0] = -gamma**i
		Bmat[i,i] = 1
	return Bmat.LLL()

def pmns(psize):
	p = random_prime(2**psize)
	# p = 2**255-19
	print("\np =",p)
	phi = 2**64
	n = floor(log(p,2)/64) + 1
	print("\nn de départ =",n)
	lamb_max = 200
	
	while True:
		lamb = 2

		while 2*n*abs(lamb)*p**(1/n) < phi and lamb < lamb_max:
			k = GF(p)
			pol = PolynomialRing(k,"X")
			X = pol("X")
			E = X**n - lamb
			E = ZZ["X"](list(E))
			Xp = sqandmult(X,p,E)
			d = xgcd((Xp-X)%E,E)[0]

			if d == 1:
				lamb += 2
				continue

			gamma = int(-d[0]%p)
			E = ZZ["X"](f"X^{n}-{lamb}")

			Bmat = matrice_base_courte(p,n,gamma)

			Bnorm = Bmat.norm(1)
			if 2*n*abs(lamb)*Bnorm < phi:
				rho = int(Bnorm - 1)
			else :
				lamb += 2
				continue
			# print("\nbase courte :\n",Bmat)
			# print("\nrho =", rho, "\ngamma =", gamma)

			liste = None
			for ligne in Bmat :
				if ligne[0]%2 == 1:
					liste = ligne
					break
			# print("\nListe retenue :",liste)

			M = ZZ["X"](list(liste))
			# print("\nM(X) =",M)

			dM,uM,_ = xgcd(M,E)
			dM = int(dM)
			dinv = pow(dM,-1,phi)
			Minv = dinv*uM%phi
			Minv = ZZ["X"](Minv)
			# print("\nM^(-1)(X) =",Minv)

			Acoeffs = [random.randrange(-rho+1,rho) for _ in range(n)]
			Bcoeffs = [random.randrange(-rho+1,rho) for _ in range(n)]
			A = ZZ["X"](Acoeffs)
			B = ZZ["X"](Bcoeffs)
			# print("\nA(X) =",A,"\nB(X) =",B)

			C = (A*B)%E
			C = ZZ["X"](C)
			# print("\nC(X) =",C)

			Q = ((C*Minv)%E)%phi
			Q = ZZ["X"](Q)
			# print("\nQ(X) =",Q)

			Cprime = (C - (Q*M)%E)/phi
			Cprime = ZZ["X"](Cprime)
			# print("\nC'(X) =",Cprime)

			print("\nn =", n, "\nlambda =", lamb)

			print("\nEst-ce que ça fonctionne :", Cprime(gamma)%p==(A(gamma)*B(gamma)*pow(phi,-1,p))%p,"\n")

			return
		n = n+1

pmns(2048)
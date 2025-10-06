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

def pmns():
	p = 2**255 - 19
	n = 5
	lambd = 2
	phi = 2**64
	K = GF(p)
	pol = PolynomialRing(K, "X")
	X = pol("X")
	E = X**n - 2
	E = ZZ["X"](f"X^{n}-{lambd}")
	Xp = sqandmult(X,p,E)
	d, u, v = xgcd((Xp-X)%E, E)
	gamma = -d[0]
	tab = [[-(gamma)**j]+[1 if i==j else 0 for i in range(1,n)] for j in range(1,n)]
	tab = [[p] + [0 for _ in range(n-1)]] + tab
	B = matrix(ZZ, tab).LLL()
	# print(B)
	if 2*n*abs(lambd)*B.norm(1) < phi :
		rho = int(B.norm(1) - 1)
	print(p)

	M = 0
	for T in B :
		M = T if T[0]%2 else M
		break
	

	M = ZZ["X"](list(M))

	dM, uM, vM = xgcd(M, E)
	dM = int(dM)
	dm1 = pow(dM, -1, phi)

	Mm1 = dm1 * uM%phi
	Mm1 = ZZ["X"](Mm1)
	print(f"M-1(X)={Mm1}")

	coeffs_A = [randrange(-rho+1,rho) for _ in range(n)]
	A = ZZ["X"](coeffs_A)
	print(f"A(X)={A}")
	coeffs_B = [randrange(-rho+1,rho) for _ in range(n)]
	B = ZZ["X"](coeffs_B)
	print(f"B(X)={B}")

	gamma = int(gamma)

	C = (A*B)%E
	print(f"C(X)={C}")

	Q = ((C*Mm1)%E)%phi
	print(f"Q(X)={Q}")

	Cp = (C - (Q*M)%E)/phi
	Cp = ZZ["X"](Cp)
	print(f"C'(X)={Cp}")

	result = (Cp(gamma)%p == (A(gamma)*B(gamma)*pow(phi,-1,p))%p)
	print(f"Result : {result}")

pmns()
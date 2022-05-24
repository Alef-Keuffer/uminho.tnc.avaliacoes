# By Alef Keuffer

# aplicar o CRT: use `crt` function in SageMath
# Euler for calculating aᵉ (mod n): @power_mod(a,e,n) see [Stein,pp.34-36 §2.3.2]

# find primitive root of natural: @primitive_root(n)


# Find last n digits of aᵉ in base b by hand ... calculate aᵉ (mod bⁿ).
# 1. a^e' ← a^e (mod b^n)
# 2. ...

def crack_when_pq_close(n):
    """When p and q are Close (uses Fermat Factorization Method) p.62"""
    t = Integer(ceil(sqrt(n)))
    while True:
        k = t^2 - n
        if k > 0:
            s = Integer(int(round(sqrt(t^2 - n))))
            if s^2 + n == t^2:
                return t+s, t-s
        t += 1

def pollard_rho(n,g=lambda x: mod(x**2+1,n)):
    """Cormen, Thomas H.; Leiserson, Charles E.; Rivest, Ronald L. & Stein, Clifford (2009). "Section 31.9: Integer factorization". Introduction to Algorithms (third ed.). Cambridge, MA: MIT Press. pp. 975–980. ISBN 978-0-262-03384-8. (this section discusses only Pollard's rho algorithm). Bu taken from wikipedia.

    It uses only a small amount of space, and its expected running time is proportional to the square root of the size of the smallest prime factor of the composite number being factorized.

    The algorithm takes as its inputs n, the integer to be factored; and g(x), a polynomial in x computed modulo n. In the original algorithm,

        g(x)=x^2-1  (mod n),

    but nowadays it is more common to use

        g(x)=x^2+1  (mod n)

    """

    x = 2
    y = 2
    d = 1

    while d == 1:
        x = g(x)
        y = g(g(y))
        d = gcd(x - y, n)

    if d == n:
        return None
    else:
        return d

def pollard_p_minus_one(N, B=10^5, stop=10):
    """Pollard’s (p − 1)-Method

    Given a positive integer N and a bound B, this algorithm attempts to find a
    nontrivial factor g of N . (Each prime p | g is likely to have the property
    that p − 1 is B-power smooth.)
    """
    m = prod([p^int(math.log(B)/math.log(p))
              for p in prime_range(B+1)])
    for a in [2..stop]:
        x = (Mod(a,N)^m - 1).lift()
        if x == 0: continue
        g = gcd(x, N)
        if g != 1 or g != N: return g
    return 1

def crack_rsa_given_phi_n(n, phi_n):
    """Factoring n Given φ(n)"""
    R.<x> = PolynomialRing(QQ)
    f = x^2 - (n+1 - phi_n)*x + n
    return [b for b, _ in f.roots()]

# cifrar e decifrar com o RSA

def rsa(bits):
    # from W. A. Stein, Elementary number theory: primes, congruences, and secrets a computational approach. New York, NY: Springer, 2009. p. 58.
    # only prove correctness up to 1024 bits
    proof = (bits <= 1024)
    p = next_prime(ZZ.random_element(2**(bits//2+1)),
            proof=proof)
    q = next_prime(ZZ.random_element(2**(bits//2+1)),
            proof=proof)
    n = p*q
    phi_n = (p-1)*(q-1)
    while True:
        e = ZZ.random_element(1,phi_n)
        if gcd(e,phi_n) == 1: break
    d = lift(Mod(e,phi_n)^(-1))
    return d, n, e

def rsa_e(bits,e):
    # only prove correctness up to 1024 bits
    proof = (bits <= 1024)
    while True:
        p = next_prime(ZZ.random_element(2**(bits//2+1)),
                proof=proof)
        q = next_prime(ZZ.random_element(2**(bits//2+1)),
                proof=proof)
        n = p*q
        phi_n = (p-1)*(q-1)
        if gcd(e,phi_n) == 1: break
    d = lift(Mod(e,phi_n)^(-1))
    return d, n, e

def rsa_encrypt(m, n, e):
    assert m < n # message must be in ℤ/nℤ
    return lift(Mod(m,n)^e)

def rsa_decrypt(c, d, n):
    return lift(Mod(c,n)^d)

# encontrar/verificar uma raíz primitiva de um natural

def multiplicative_group(n):
    return {IntegerModRing(n)(e)
            for e in IntegerModRing(n).list_of_elements_of_multiplicative_group()}

def findOrders(n):
    return {(e,e.multiplicative_order())
            for e in multiplicative_group(n)}

def all_prim_roots(n):
    return {e[0] for e in findOrders(n) if e[1] == euler_phi(n)}

def num_prim_roots(n):
    return euler_phi(euler_phi(n))

def is_prim_root(a,n):
    return a in all_prim_roots(n)

# By Alef Keuffer

# Mini-teste 1

# Considere
n=1431702961131339621602945101050088245802771644687686106472404992618656781417216251296575852739319584928520070989102317406096771341195807540255351846662720304620283140579998725212022478694288443,
# o produto de dois primos p e q, com p<q, é
phi_n=1431702961131339621602945101050088245793639679386713169093313712250573769060652354811397118799813658291947433707988011985679885048571032964071342275971679191493672348651555118274559858896490208
# entao p,q sao
crack_rsa_given_phi_n(n,phi_n)
---
# De uma chave pública RSA
(n,e)=(5746373959926501631, 816210188025244429)
# sabe-se que
p=1434802597
# divide n. Então o expoente de decifração, d, é igual a
q = n/p
phi_n = (p-1)*(q-1)
d = inverse_mod(e,phi_n)
---
# Considere o natural n=121
7 in all_prim_roots(121)
---
# Para
n=7*11
# e
a=60*10^6+3
# , o resto da divisão inteira de 3^a por n é
Mod(3^a,n)
---
# Dada a chave pública RSA com
(n,e) = (1020503581,180933103)
# , a cifração de
mens=1234
# é
rsa_encrypt(mens,n,e)
---
# Considere
n=2^30-1
m=2^29-1
# e o isomorfismo f de anéis Zₙₘ para Zₙ x Zₘ. A imagem recíproca de
(a,b) = (402374628, 220354304)
# é
crt([a,b],[n,m])
# coding=utf-8
# By Alef Keuffer

# cifrar/decifrar usando o ElGamal
def elgamal_encrypt(public_key, message, k_random=None):
    # b_power is the result of an exponentiation operation (r^α where α is a secret)
    p_prime, r_primitive_root, b_power = public_key
    Zp = IntegerModRing(p_prime)
    r_primitive_root, b_power = Zp(r_primitive_root), Zp(b_power)
    if k_random is None:
        k_random = randint(2, p_prime - 2)
    gama, delta = r_primitive_root ^ k_random, message * b_power ^ k_random
    return gama, delta


def elgamal_decrypt(public_key, a_exponent, cryptogram):
    # b_power is the result of an exponentiation operation (r^α where α is a secret)
    p_prime, r_primitive_root, b_power = public_key
    Zp = IntegerModRing(p_prime)
    gama, delta = cryptogram
    gama, delta = Zp(gama), Zp(delta)
    decrypted_message = (gama ^ a_exponent) ^ (-1) * delta
    return decrypted_message


# reconhecer raizes primitivas

def multiplicative_group(n):
    return {IntegerModRing(n)(e)
            for e in IntegerModRing(n).list_of_elements_of_multiplicative_group()}


def findOrders(n):
    return {(e, e.multiplicative_order())
            for e in multiplicative_group(n)}


def all_prim_roots(n):
    return {e[0] for e in findOrders(n) if e[1] == euler_phi(n)}


def num_prim_roots(n):
    return euler_phi(euler_phi(n))


def is_prim_root(a, n):
    return a in all_prim_roots(n)


# verificar o valor do índice (aka logaritmo discreto)
# b ≡ rᵃ (mod n) ⟺ indᵣb ≡ a (mod n)
def is_index_of(r, b, a, n):
    """
    returns true if a (mod n) is index of b in base r, false otherwise
    """
    return power_mod(r, a, n) == b


def index_of(n, r, b):
    """
    returns indᵣb
    """
    return discrete_log(b, IntegerModRing(n)(r))


# aplicar o Lema de Gauss para calcular o símbolo de Legendre

"""
Let p=17 and a=7.

There are 16 nonzero elements [1...16].

Consider the first half [1...8] and multiply them all by 7 to get

    7, 14, 4, 11, 1, 8, 15, 5.

We single out 14, 11 and 15 because they are greater than p/2 (that is, 9 or higher).

So exactly 3 of them are greater than p/2.

Gauss' Lemma states that if we take this 3 and raise -1 to this power, then we have legendre_symbol(a,b)

that is:

    legendre_symbol(7,17) = (-1)³ = -1
"""

# Identificar pseudo-primos de Euler
# (ou seja, que consigam afirmar se um certo número apresentado é o psprimo Euler de base dada).
# Em particular que não se esquecessem que, tal como todos os pseudo-primos,
# os de Euler são necessariamente números compostos.

"""
an odd composite integer n is called an Euler pseudoprime to base a, if a and n are coprime, and

a^((n-1)/2) ≡ ± 1 (mod n)
"""


# sagemath is_pseudoprime() may be the same thing

def is_euler_pseudoprime(n, base=2):
    assert n % 2 == 1
    return gcd(a, n) == 1 and (power_mod(a, (n - 1) // 2, n) == IntegerModRing(1) or power_mod(a, (n - 1) // 2,
                                                                                               n) == -1) and not is_prime(
        n)


# Perguntas do Mini-teste 2

# Pergunta 1
# Dada a chave pública ElGamal
pub_key = (59879295262580794019, 2, 46532948489070777896)  # (p,r,b)
# sabendo que o índice de 46532948489070777896 na base 2 é igual a
a = 26889797028840904448
# a decifração de
cryptogram = (31508193067819085597, 42769957659645449029)
# é igual a ___.
elgamal_decrypt(pub_key, a, cryptogram) #1342

# Pergunta 2
# Considere o primo
p = 17
# e o natural
a = 5
# Então o número dos menores resíduos de ka maiores que p/2, com k entre 1 e (p-1)/2 é igual a ___.
len([Mod(k * a, p) for k in range(1, (p - 1) / 2 + 1) if Mod(k * a, p) > p // 2])  # == len([10,15,13]) == 3

# Pergunta 3
# Considere o primo
p = 16456780161579311951
# e a raiz primitiva
r = 10632721783361421121
# de p. O índice de 12124670967023022731 na base r é igual a
[a for a in
 [1945473589762547835, 9108492322933321998, 765875512, 78664311121]
 if is_index_of(r, 12124670967023022731, a, p)]  # == [9108492322933321998]

# Pergunta 4
# São pseudo primos de Euler de base 5
[p for p in [99, 781, 101, 1541] if is_euler_pseudoprime(p, 5)]  # == [781, 1541]

# Pergunta 5
# Indique todos os elementos que são raiz primitiva de 233
[r for r in [21, 86, 46, 207] if r in all_prim_roots(233)]  # == [21, 86]

# Pergunta 6
# Dada a chave pública ElGamal
pub_key = (39957963614327378639, 13, 2518634842003570977)
# a cifração de
mens = 1234
# com o parâmetro aleatório
k = 4321
# é
elgamal_encrypt(pub_key, mens, k)  # == (11727766830592979831, 35640834212194464972)

import random
from time import time
from statistics import mean

print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP9 : FACTORISATION DES ENTIERS                                             #
# *************************************************************************** #
# *************************************************************************** #
""")

# CONSIGNES
#
# Les seules lignes a modifier sont annoncee par "Code pour l'exercice"
# indique en commmentaire et son signalees
# Ne changez pas le nom des variables
#
# CONSEILS
#
# Ce modele vous sert a restituer votre travail. Il est deconseille d'ecrire
# une longue suite d'instruction et de debugger ensuite. Il vaut mieux tester
# le code que vous produisez ligne apres ligne, afficher les resultats et
# controler que les objets que vous definissez sont bien ceux que vous attendez.
#
# Vous devez verifier votre code en le testant, y compris par des exemples que
# vous aurez fabrique vous-meme.
#


#reset()
print("""\
# ****************************************************************************
# DIVISEURS SUCCESSIFS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=2*3*3*5*5*5*7*11*11

# Code pour l'EXERCICE

def div_successives(n):
    f = []
    d = 2
    while d^2<=n:
        cnt_d = 0
        while n.mod(d)==0:
            n = n//d
            cnt_d+=1
        if cnt_d>0:
            f.append((d,cnt_d))
        d+=1
    if n!=1:
        f.append((n,1))
    return f

# # Affichage des resultats

print(div_successives(n))
for n in range(2,10):
    assert(div_successives(ZZ(n))==list(factor(ZZ(n))))


#reset()
print("""\
# ****************************************************************************
# FACTORISATION D'UN NOMBRE B-FRIABLE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=2*3*3*5*5*5*7*11*11
P=[p for p in primes(12)]

# Code pour l'EXERCICE

def div_successives_friable(n, P):
    f = []
    for d in P:
        cnt_d = 0
        while n.mod(d)==0:
            n = n//d
            cnt_d+=1
        if cnt_d>0:
            f.append((d,cnt_d))
        d+=1
    if n!=1:
        f.append((n,1))
    return f

# # Affichage des resultats

print(div_successives_friable(n,P))
for n in range(2,10):
    assert(div_successives_friable(ZZ(n),P)==list(factor(ZZ(n))))


#reset()
print("""\
# ****************************************************************************
# RHO DE POLLARD
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=222763

# Code pour l'EXERCICE

def myPollardrho(n,max_iter=100):
    x0=Integers(n).random_element()
    y0 = x0
    Poln.<x> = PolynomialRing(Integers(n))
    f = x^2+1
    x0 = f(x0)
    y0 = f(f(y0))
    g = (x0-y0).gcd(n)
    i = 1
    while g<=1 and i<max_iter:
        x0 = f(x0)
        y0 = f(f(y0))
        g = (x0-y0).gcd(n)
        i+=1
    if g.mod(n) in [0,1]:
        return 1
    else:
        return g



# # Affichage des resultats

print(myPollardrho(n))

for _ in range(5):
    n=ZZ.random_element(3,100)
    print(n, 
      "| Resultat rho de Pollard : ", 
      myPollardrho(n), 
      " | n est-il composé ?",not n.is_prime())




#reset()
print("""\
# ****************************************************************************
# P-1 DE POLLARD
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=1323269
b = 11

# Code pour l'EXERCICE

# def myPollardpm1(n,b=15):
#     a = Integers(n).random_element()
#     g = a.gcd(n)
#     if g>1:
#         return g
#     B = lcm(range(1,b+1))
#     for p in prime_range(2,B):
#         p = Integer(p)
#         if p.is_prime():
#             alpha = 0
#             while p^(alpha+1)<=B:
#                 alpha+=1
#             #for _ in range(1,alpha+1):
#             #   a = a^p.mod(n)
#             a = a^(p^alpha).mod(n)
#     g = (a-1).gcd(n)
#     if 1<g and g<n:
#         return g
#     else:
#         return 1


def myPollardpm1(n,b=50):
    #a = Integers(n).random_element()
    a = random.randrange(2,n)
    g = gcd(a,n)
    if g>1:
        return g

    for p in prime_range(2,b+1):
        alpha = 1
        while p**(alpha+1)<=b:
           alpha+=1
        #a = a^(p^alpha)

        a = pow(a,p**alpha,n)
        g = gcd(a-1, n)
        if 1<g<n:
            return g

    return 1



# # Affichage des resultats

print(myPollardpm1(n,b=50))

for _ in range(5):
    n=ZZ.random_element(3,100)
    print(n, 
      "| Resultat rho de Pollard : ", 
      myPollardpm1(n), 
      " | n est-il composé ?",not n.is_prime())



nmin, nmax = 400, 600
nbtests = 20

freqs_p_1 = []
times_p_1 = []
freqs_p = []
times_p = []
for n in range(nmin, nmax):
    if is_prime(n):
        freqs_p.append(0)
        times_p.append(0)
        freqs_p_1.append(0)
        times_p_1.append(0)
        continue
    successes_p = 0
    successes_p_1 = 0
    ts_p = []
    ts_p_1 = []
    for i in range(nbtests):
        t0 = time()
        res_p_1 = myPollardpm1(n, 50)
        ts_p_1.append(time() - t0)
        t1 = time()
        res_p = myPollardrho(n)
        ts_p.append(time() - t1)
        if res_p > 1:
            successes_p += 1
        if res_p_1 > 1:
            successes_p_1 += 1
    freqs_p.append(successes_p/nbtests)
    times_p.append(mean(ts_p))
    freqs_p_1.append(successes_p_1/nbtests)
    times_p_1.append(mean(ts_p_1))


# Plot Success Frequency
bar_chart(
    freqs_p_1,
    title='Success Frequency Pollard p-1 vs 400<=n<=600',
    figsize=[8,4]
).save('success_pollard_p-1.png')
print("saved 'success_pollard_p-1.png'")

# Plot Average Runtime
bar_chart(
    times_p_1,
    title='Average Runtime (s) Pollard p-1 vs 400<=n<=600',
    figsize=[8,4]
).save('runtime_pollard_p-1.png')
print("saved 'runtime_pollard_p-1.png'")

bar_chart(
    freqs_p,
    title='Success Frequency Pollard rho vs 400<=n<=600',
    figsize=[8,4]
).save('success_pollard_pho.png')
print("saved 'success_pollard_pho.png'")

# Plot Average Runtime
bar_chart(
    times_p,
    title='Average Runtime (s) Pollard rho vs 400<=n<=600',
    figsize=[8,4]
).save('runtime_pollard_pho.png')
print("saved 'runtime_pollard_pho.png'")





#reset()
print("""\
# ****************************************************************************
# CRIBLE QUADRATIQUE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=2886

# Code pour l'EXERCICE


def tonelli(p,a):
    if p.mod(4)==3:
        x = a^((p+1)//4).mod(p)
    elif p.mod(8)==5:
        x = a^((p+3)//8).mod(p)
        c = x^2.mod(p)
        if c.mod(p)!=a:
            x = x*(2^((p-1)//4)).mod(p)
    else:
        d = 2
        while d<p and legendre_symbol(d,p)!=-1:
            d+=1
        s = (p-1).valuation(2)
        t = (p-1)//2^s
        A = a^t.mod(p)
        D = d^t.mod(p)
        m = 0
        for i in range(s):
            if (A*D^m)^(2^(s-1-i)).mod(p)==-1:
                m+=2^i
        x = a^((t+1)//2)*D^(m//2).mod(p)
    return x

print("Test de racine carrée de Tonelli: ")
for p in prime_range(3,501):
    for a in Integers(p):
        if legendre_symbol(a,p)==1:
            assert(tonelli(p,a)^2.mod(p)==a)
print("réalisé avec succès")



def hensel(f,r,p,k):
    for i in range(k):
        x = f(r)*p^(-2^i).mod(p^(2^i))
        z = f.derivative()(r)^-1.mod(p^(2^i))
        y = -x*z.mod(p^(2^i))
        r = r+y*p^(2^i)
    return r


def cribleBasique(K,f,B):
    F = [1 for _ in range(K)]
    P = prime_range(B+1)
    maxf = max(abs(f(k)) for k in range(1, K+1))
    for p in P:
        alpha = 1
        #while p^(alpha+1)<=maxf:
        #    alpha+=1
        alpha = floor(log(maxf, p))
        for a in range(1,alpha+1):
            pa = p**a
            R = IntegerModRing(pa)
            P.<x> = PolynomialRing(R)
            pho_list = P(f).roots(multiplicities=False)
            for pho in pho_list:
                for k in range(pho, K+1, pa):
                    F[k-1] *= p
    cand = []
    for k in range(K):
        if f(k+1)==F[k]:
            cand.append(k+1)
    return cand

Pol.<x> = PolynomialRing(ZZ)
f = (x+54)^2-2886
print("Crible basique de f = ",f,"11-friable :",cribleBasique(15,f,11))


def cribleBasique_recherche(K,n,B):
    Pol.<x> = PolynomialRing(ZZ)
    f = x^2-n
    return [z for z in cribleBasique(K,f,B) if z>=sqrt(n)]


def cribleQuadratique(n):
    B = ceil(math.exp(1/2*(sqrt(log(n)*log(log(n))))))
    F = prime_range(2,B+1)
    m = len(F)
    S = []
    cnt = 0
    l = cribleBasique_recherche(200,n,B)
    for x in l:
        if cnt>=m+1:
            break
        S.append((x,x**2-n))
        cnt+=1
    S_fact = []
    for (x,a) in S:
        S_fact.append((x,div_successives_friable(a,F)))
    r = len(S_fact)
    M = Matrix(GF(2), m, r)
    for j, p in enumerate(F):
        for i, (x, fac) in enumerate(S_fact):
            e = 0
            for (q,exp) in fac:
                if q==p:
                    e = exp
                    break
            M[j,i] = e%2
    kernel_vecs = M.right_kernel().basis()
    if kernel_vecs==[]:
        return 1
    v = next(v for v in kernel_vecs if not v.is_zero())
    K = [i for i, bit in enumerate(v) if bit == 1]
    Z = prod(S_fact[i][0] for i in K) % n
    Y = 1
    for j, p in enumerate(F):
        tot = sum(
            exp
            for i in K
            for (q, exp) in S_fact[i][1]
            if q == p
        )
        Y = (Y*power_mod(p, tot//2,n))%n
    d = gcd(Z-Y,n)
    if d in (1, n):
        d = gcd(Z+Y,n)
    return d


# # Affichage des resultats

print("Crible quadratique pour n=2886:",cribleQuadratique(3652))

nmin, nmax = 400, 500
nbtests = 2

freqs_qs = []
times_qs = []

for n in range(nmin, nmax):
    if is_prime(n):
        freqs_qs.append(0)
        times_qs.append(0)
        continue

    succ_qs = 0
    ts_qs   = []
    for _ in range(nbtests):
        t0 = time()
        d  = cribleQuadratique(n)
        dt = time() - t0
        ts_qs.append(dt)
        if 1<d<n:
            succ_qs += 1

    freqs_qs.append(succ_qs/nbtests)
    times_qs.append(mean(ts_qs))

bar_chart(
    freqs_qs,
    title='Success Frequency Crible Quadratique vs 400<=n<=500',
    figsize=[8,4]
).save('success_cq.png')
print("Saved 'success_cq.png'")

bar_chart(
    times_qs,
    title='Average Runtime (s) Crible Quadratique vs 400<=n<=500',
    figsize=[8,4]
).save('runtime_cq.png')
print("Saved 'runtime_cq.png'")



from tp4_daboust import *

import random
print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP5 : FACTORISATION COMPLETE DE POLYNOMES UNIVARIEES                        #
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
# BERLEKAMP
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

F3 = FiniteField(3)
Pol3.<x> = PolynomialRing(F3)
f = x^3 - x^2 - 1

# Code pour l'EXERCICE

x3 = list(Pol3(x^3).mod(f))
x6 = list(Pol3(x^6).mod(f))+[0]
Q = matrix(F3,3,[[1,0,0],x3,x6]).transpose()

b1 = vector(F3,[1,0,0])
b2 = vector(F3,[0,1,1])

L = (Q-matrix.identity(3)).transpose().right_kernel().basis()

b1 = vector(F3,L[0])
b2 = vector(F3,L[1])

f0 = (Pol3(list(b2))+0).gcd(f)
f1 = (Pol3(list(b2))+1).gcd(f)
f2 = (Pol3(list(b2))+2).gcd(f)

def myB(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [f]
    Fq = Pol.base_ring()
    n = f.degree()

    Q = matrix(Fq,n)
    for j in range(n):
        elt = [Fq(c) for c in Pol(x^(j*q)).mod(f)]

        nb_zero = n-len(elt)
        elt = elt + nb_zero*[0]
        for i in range(n):
            Q[i,j] = Fq(elt[i])

    K = (Q-matrix.identity(n)).right_kernel().basis()
    K_list = []
    for i in range(len(K)):
        K_list.append(vector(Fq,K[i]))

    F = set()
    F.add(f.monic()) # polynôme unitaire (on multiplie par l'inverse du coeff dominant)
    j = 0
    while len(F)<len(K_list):
        F = set(f.monic() for f in F) # pour polynôme unitaire
        j+=1
        C = [f_t for f_t in F if f_t.degree()>1]
        for f_t in C:
            B = set()
            for e in Fq:
                a = f_t.gcd(Pol(list(K_list[j]))-e)
                if a.degree()>=1:
                    B.add(a.change_ring(Fq))
            F.remove(f_t)
            F = F.union(B)

    assert(Set(F) == Set(g for g,_ in list(f.factor())))
    return F

def random_polynomial_sans_carre(max_degree=5, max_field_size=5, timeout=1000):
    p = random.choice(list(prime_range(2, max_field_size)))
    n = random.randint(1, 4)
    q = p**n
    Fq = FiniteField(q)
    R = PolynomialRing(Fq, 'x')
    x = R.gen()
    random_poly = R.random_element(degree=max_degree)
    step=1
    while random_poly.gcd(random_poly.derivative())!=1 and step<timeout:
        random_poly = R.random_element(degree=max_degree)
        step+=1
    return random_poly


sample = [random_polynomial_sans_carre() for _ in range(100)] # 100 polynômes de corps finis sans facteur carré

test1 = [Set(myB(f)) for f in sample] == [Set(g for g,_ in list(f.factor())) for f in sample]


def myFactor(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
    f_sans_carre = myFsFC(f)
    for (g,e) in f_sans_carre:
        for h in myB(g):
            retour.append((h,e))
    assert(Set(retour) == Set(list(f.factor())))
    return retour


def random_polynomial(max_degree=10, max_field_size=10):
    p = random.choice(list(prime_range(2, max_field_size)))
    n = random.randint(1, 3)
    q = p**n
    Fq = FiniteField(q)
    R = PolynomialRing(Fq, 'x')
    x = R.gen()
    random_poly = R.random_element(degree=max_degree)
    return random_poly


sample = [random_polynomial() for _ in range(100)] # 100 polynômes de corps finis random
test2 = [Set(myFactor(f)) for f in sample] == [Set(list(f.factor())) for f in sample]


# # Affichage des resultats


print("\n$1a/ x^3 vaut",Pol3(x3)," et x^6 vaut",Pol3(x6))
print("La matrice de Petr Berlekamp est")
print(Q)

print("\n$1b/ On a Q * b1 - b1 = ")
print(Q*b1-b1)
print("et Q * b2 - b2 = ")
print(Q*b2-b2)
print("avec b1 = ",b1,"et b2 = ",b2)
print("Les facteurs obtenus pour f sont ",myB(f))
print("Test de myB sur 100 examples :", test1)

print("Factorisation de f :", myFactor(f))
print("Test de factorisation sur 100 examples :", test2)

#reset()
print("""\
# ****************************************************************************
# RELEVEMENT DE HENSEL
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

PolZZ.<x> = PolynomialRing(ZZ)
m = 5
f = x^4-1
g = x^3+2*x^2-x-2
h = x-2
d,ss,tt = xgcd(g,h)
s=PolZZ(ss/mod(d,m)); t=PolZZ(tt/mod(d,m))

# m = 3
# f = x^4+78*x^3-2556*x^2+4389*x-722
# g = x^2+x-1
# h = x^2-x-1
# s = -x+1
# t = x+1

# Code pour l'EXERCICE



def polynomeCentre(f):
    Pol=f.parent()
    x=Pol.gen()
    retour = [coeff.lift_centered() for coeff in f.list()]
    return PolZZ(retour) # permet de garder le centrage et pas rebasculer en forme canonique


def myHensel(f,g,h,s,t,m):
    Pol=f.parent()
    x=Pol.gen()
    Zm = Integers(m)
    Polz = PolynomialRing(Zm,'x')
    Zm2 = Integers(m^2)
    Polz2 = PolynomialRing(Zm2,'x')
    prod_gh = Polz(g*h)
    assert(Polz(f-prod_gh)==0)
    prod_euclide = polynomeCentre(Polz(s*g+t*h))
    assert(prod_euclide==Polz(1))
    assert(s.degree()<h.degree() and t.degree()<g.degree())

    f = Polz2(f)
    g = Polz2(g)
    h = Polz2(h)
    s = Polz2(s)
    t = Polz2(t)

    e = Polz2(f-g*h)
    q,r = (s*e).quo_rem(h)
    g_star = polynomeCentre(Polz2(g+t*e+q*g))
    h_star = polynomeCentre(Polz2(h+r))
    b = s*g_star+t*h_star-1
    c,d = (s*b).quo_rem(h_star)
    s_star = polynomeCentre(Polz2(s-d))
    t_star = polynomeCentre(Polz2(t-t*b-c*g_star))

    retour = g_star,h_star,s_star,t_star
    return retour


def myHenselItere(f,g,h,s,t,m,l):
    Pol=f.parent()
    x=Pol.gen()
    Zm = Integers(m)
    Polz = PolynomialRing(Zm,'x')
    Zm2 = Integers(m^2)
    Polz2 = PolynomialRing(Zm2,'x')
    prod_gh = polynomeCentre(Polz(g*h))
    assert(g.gcd(h)==Pol(1))
    assert(Polz(f-prod_gh)==0)

    k = 0
    while 2^k<l:
        g,h,s,t = myHensel(f,g,h,s,t,m^(2^k))
        k+=1
    retour = g,h,s,t
    return retour

reponseQ5="\nSi f se décompose en n facteurs, appliquer récursivement le relèvement aux (n-1) facteurs g, et h le dernier facteur, en ."

# # Affichage des resultats

print("\n$1b/ Relèvement de ",f,"= (",g,")*(",h,")")
print(myHensel(f,g,h,s,t,m))
print("\nSur plusieurs étapes de ",f," :",myHenselItere(f,g,h,s,t,m,6))
print(reponseQ5)


#reset()
print("""\
# ****************************************************************************
# FACTORISATION AVEC LLL
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p=13
k=4
m=p^k

PolZZ.<x> = PolynomialRing(ZZ)
f = x^4-x^3-5*x^2+12*x-6
Fp = Integers(p)
PolFp.<x> = PolynomialRing(Fp)
Fp2 = Integers(p^2)
PolFp2.<x> = PolynomialRing(Fp2)
f_p = PolFp(f)
f_p_fact = myFactor(f_p) # on a 4 facteurs de degré 1

alpha=-f_p_fact[0][0](0)
beta=-f_p_fact[1][0](0)
gamma=-f_p_fact[2][0](0)
delta=-f_p_fact[3][0](0)

g1 = PolFp(f_p_fact[0][0]*f_p_fact[1][0]) # on split f en deux facteurs de degré 2
h1 = PolFp(f_p_fact[2][0]*f_p_fact[3][0])
d,ss,tt = xgcd(g1,h1)
s1 = PolFp(ss//mod(d,p)); t1 = PolFp(tt/mod(d,p))
g2,h2,s2,t2 = myHenselItere(f,g1,h1,s1,t1,p,4) # on a deux facteurs de degré 2, réduits modulo 13^4

g3 = PolFp(f_p_fact[0][0]) # on a g2 qui est réduit modulo 13 aux deux premiers des quatre facteurs initiaux de f: on le relève modulo 13^4
h3 = PolFp(f_p_fact[1][0])
d,ss,tt = xgcd(g3,h3)
s3 = PolFp(ss//mod(d,p)); t3 = PolFp(tt/mod(d,p))
g3,h3,s3,t3 = myHenselItere(g2,g3,h3,s3,t3,p,4) # on a les deux premiers facteurs de f réduits modulo 13^4

g4 = PolFp(f_p_fact[2][0]) # pareil pour h2
h4 = PolFp(f_p_fact[3][0])
d,ss,tt = xgcd(g4,h4)
s4 = PolFp(ss//mod(d,p)); t4 = PolFp(tt/mod(d,p))
g5,h5,s5,t5 = myHenselItere(h2,g4,h4,s4,t4,p,4) # les deux derniers réduits modulo 13^4


alphahat=-g3(0)
betahat=-h3(0)
gammahat=-g5(0)
deltahat=-h5(0)

print("\nVérification racines en réduction modulo 13^4 ne le sont pas sur Z[x] :",f(alphahat),"!=0",f(betahat),"!=0",f(gammahat),"!=0",f(deltahat),"!=0")

PolZZ.<x> = PolynomialRing(ZZ)
Fm = Integers(m)
PolFm.<x> = PolynomialRing(Fm)

u1 = x+7626
print("\nChoisissons une première racine -u1(0) =",-u1(0),"telle que f s'annule bien f(-u1(0)) =",PolFm(f(-u1(0))),"mod 13⁴")
j = 3
d = u1.degree()

u2 = x+9864
print("\nChoisissons la deuxième racine -u2(0) =",-u2(0),"telle que f s'annule bien f(-u2(0)) =",PolFm(f(-u2(0))),"mod 13⁴")

# Code pour l'EXERCICE

def find_fact(u):
    p0 = list(u)+[0]*(j-len(list(u)))
    p1 = u*x
    p1 = list(p1)
    p1 = p1 + [0]*(j-len(p1)) # u*x^i
    M = matrix(ZZ,j,j+(j-d))

    for i in range(j):
        M[i,0] = p0[i]
        M[i,1] = p1[i]

        for k in range(j-d,2*j-d):
            for i in range(j):
                if i==(k-j+d):
                    M[i,k] = m
                else:
                    M[i,k] = 0

    L = (M.transpose().LLL().transpose())
    l = [L.column(j) for j in range(L.ncols()) if vector(L.column(j))!=vector([0,0,0])]
    g = PolZZ(list(l[0])) # court vecteur

    return l,g,PolZZ(f).gcd(g)

# # Affichage des resultats

print("\n$1a/ Les racines sont", alpha, beta, gamma, delta,"modulo",p)
print("\n$1b/ Les racines sont", alphahat, betahat, gammahat, deltahat,"modulo",m)

fact1 = find_fact(u1)
print("\nPremier facteur")
print("Base LLL-réduite du réseau :",fact1[0])
print("Vecteur court de la base LLL-réduite :",fact1[1])
print("Un facteur de f est alors gcd(f,g1) =",fact1[2])

fact2 = find_fact(u2)
print("\nDeuxième facteur")
print("Base LLL-réduite du réseau :",fact2[0])
print("Vecteur court de la base LLL-réduite :",fact2[1])
print("Un facteur de f est alors gcd(f,g1) =",fact2[2])

print("\nFactorisation complète de f =",f)
print(f,"=","(",fact1[2],")","(",fact2[2],")")

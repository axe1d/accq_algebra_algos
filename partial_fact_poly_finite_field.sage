import random

print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP4 : FACTORISATION DE POLYNOMES UNIVARIEES SUR CORPS FINIS                 #
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
# FACTORISATION DES PUISSANCES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

F1849.<omega> = FiniteField(43^2,modulus=x^2+1)
Pol1849.<x> = PolynomialRing(F1849)
f=x^172+(3-2*omega)*x^129-5*omega*x^86+(2 + 4*omega)*x^43-1-omega 

F9.<alpha> = FiniteField(9)
Pol9.<y> = PolynomialRing(F9)
g = y^30-y^15+alpha*y^3+1

# Code pour l'EXERCICE

def racine_p_polynome(f):
    assert(f.derivative()==0)
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    u=Pol(0)
    deg = 0
    for coeff in list(f)[0:len(list(f)):p]:
        new_coeff = coeff^(q/p)
        u = u+(new_coeff)*x^deg
        deg += 1
    assert(u^p==f)
    return u

sample = [Pol1849.random_element()(x^(Pol1849.base_ring().characteristic())) for _ in range(100)] # 100 polynômes de dérivée nulle
test = [racine_p_polynome(f)^(Pol1849.base_ring().characteristic()) for f in sample] == sample

# # Affichage des resultats

print( "\n$ Question 3")
print( "La racine de",f,"est",racine_p_polynome(f))
print( "La racine de",g,"est",racine_p_polynome(g))
print( "\n$ Question 4")
print( "Test sur 100 exemples : ",test)



#reset()
print("""\
# ****************************************************************************
# FACTORISATION SANS FACTEURS CARRES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

F7 = FiniteField(7)
Pol7.<x> = PolynomialRing(F7)
f = x^10 +6*x^9 +3*x^7 +3*x^3 +4*x^2 +2

# Code pour l'EXERCICE

def myFsFC(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [f,1]
    if f.degree()<=0:
        return []
    elif f.derivative()!=0:
        i = 1
        retour = []
        t = f.gcd(f.derivative())
        u = f//t
        while u!=1:
            y = t.gcd(u)
            if i%p!=0 and u//y!=1:
                retour.append((u//y,i))
            i+=1
            u = y
            t = t//y
        if t!=1:
            u_tilde = racine_p_polynome(t)
            retour = retour + [(s,p*i) for (s,i) in myFsFC(u_tilde)]
    else:
        u_tilde = racine_p_polynome(f)
        retour = [(s,p*i) for (s,i) in myFsFC(u_tilde)]

    assert(prod([f^e for (f,e) in retour ]) == f)
    return retour


def random_polynomial(max_degree=50, max_field_size=100):
    p = random.choice(list(prime_range(2, max_field_size)))
    n = random.randint(1, 3)
    q = p**n
    Fq = GF(q, 'a')
    R = PolynomialRing(Fq, 'x')
    x = R.gen()
    random_poly = R.random_element(degree=max_degree)
    return random_poly


sample = [random_polynomial() for _ in range(1000)] # 1000 polynômes de corps finis random
test = [prod([f^e for (f,e) in myFsFC(f)]) for f in sample] == sample

# # Affichage des resultats

print( "\n$ Question 2")
print( "La factorisation de",f,"est",myFsFC(f), "et celle de SageMath est", factor(f))
print( "\n$ Question 4")
print( "Test sur 1000 exemples : ",test)


#reset()
print("""\
# ****************************************************************************
# FACTORISATION ETAGEE EN DEGRES DISTINCTS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

F5 = FiniteField(5)
Pol.<x>=PolynomialRing(F5)
f = x^10-2*x^9+x^8+x^7-x^6-2*x^5+2*x^4+2*x^3-x

# Code pour l'EXERCICE
def exp_fast(h,n,f):
    Pol=h.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    n_bin = n.bits()
    result = Pol(1)
    base = h.mod(f)
    for bit in reversed(n_bin):
        result = (result*result).mod(f)
        if bit==1:
            result = (result*base).mod(f)
    return result


def myFEDD(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
    hi = x
    fi = f
    i = 0
    while fi!=1:
        i+=1
        hi = (exp_fast(hi,q,f)).mod(f)
        gi = (hi-x).gcd(fi)
        fi = fi//gi
        retour.append(gi)

    assert(prod(retour) == f)
    return retour


def random_polynomial_sans_carre(max_degree=10, max_field_size=10, timeout=1000):
    p = random.choice(list(prime_range(2, max_field_size)))
    n = random.randint(1, 3)
    q = p**n
    Fq = GF(q, 'a')
    R = PolynomialRing(Fq, 'x')
    x = R.gen()
    random_poly = R.random_element(degree=max_degree)
    step=1
    while random_poly.gcd(random_poly.derivative())!=1 and step<timeout:
        random_poly = R.random_element(degree=max_degree)
        step+=1
    return random_poly


sample = [random_polynomial_sans_carre() for _ in range(100)] # 100 polynômes de corps finis sans facteur carré
test = [prod(myFEDD(f)) for f in sample] == sample
#test = False

# # Affichage des resultats

print( "\n$ Question 1")
print( "La factorisation de",f,"est",myFEDD(f))
print( "Test sur 100 exemples : ",test)

#reset()
print("""\
# ****************************************************************************
# CANTOR-ZASSENHAUSS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

q=9
d=4
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq)
l = [f for f in Polq.polynomials(of_degree=d) if f.is_irreducible()
 and f.leading_coefficient()==1]

# Code pour l'EXERCICE

def choose_rd_poly(Pol, Fq, bound):
    f = []
    for _ in range(bound):
        x = Fq.random_element()
        f.append(x)
    return Polq(f)


def myCZ(f,d):
    MAX_STEPS=100
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    Fq = FiniteField(q)
    retour = []
    step = 0
    bound = 2*d
    if f.degree()<=d:
        return [f]

    u = choose_rd_poly(Pol,Fq,bound)

    b = f.gcd(exp_fast(u,(q^d-1)//2,f)-1)
    while step<MAX_STEPS:
        step+=1
        u = choose_rd_poly(Pol,Fq,bound)

        b = f.gcd(exp_fast(u,(q^d-1)//2,f)-1)
        if 0<b.degree()<f.degree():
            return myCZ(b,d) + myCZ(f//b,d)
    return [f]


def Tr(u,m,f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    h = Pol(0)
    for i in range(m):
        h = (h + exp_fast(u,2^i,f)).mod(f)
    return h


def myCZ2(f,d):
    MAX_STEPS=1000
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    bound = 2*d
    assert(q.bits().index(1)==len(q.bits())-1) # puissance de deux uniquement
    retour = []

    step = 0
    if f.degree()==d:
        return [f]

    r = random.randint(1,2*d-1)

    u = choose_rd_poly(Pol,Fq,bound)

    k = len(list(q.bits()))-1
    b = f.gcd(Tr(u,k*d,f))
    while step<MAX_STEPS:
        step+=1
        r = random.randint(1,2*d-1)
        u = choose_rd_poly(Pol,Fq,bound)
        k = len(list(q.bits()))-1
        b = f.gcd(Tr(u,k*d,f))
        if 0<b.degree()<f.degree():
            return myCZ2(b,d) + myCZ2(f//b,d)
    return [f]


def poly_random_degree_d(Pol, l, max_nb=10):
    h = Pol(1)
    n = random.randint(2,max_nb)
    n = min(n,len(l))
    chosen = random.sample(l,n)
    for hi in chosen:
        h = h*hi
    assert(h.gcd(h.derivative())==1)
    return h

#precomputed_poly = precompute_polynomials(Polq, 2*d-1)
sample2 = [poly_random_degree_d(Polq,l) for _ in range(100)]
result = [myCZ(f,d) for f in sample2]
succeed = sum([result[i]!=[sample2[i]] for i in range(len(sample2))])
test  = [prod(result[i]) for i in range(len(sample2))] == sample2

# # Affichage des resultats
print( "Test sur 100 exemples : ", test)
print("Nombre de tests terminés", succeed)

n = 2
q=2^n
d=3
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq)
l = [f for f in Polq.polynomials(of_degree=d) if f.is_irreducible()
 and f.leading_coefficient()==1]
sample3 = [poly_random_degree_d(Polq,l) for _ in range(100)]
succeed = sum([myCZ2(f,d)!=[f] for f in sample3])
test  = [prod(myCZ2(f,d)) for f in sample3] == sample3
print( "Test sur 100 exemples en caractéristique 2 : ", test)
print("Nombre de tests terminés", succeed)


#reset()
print("""\
# ****************************************************************************
# FACTORISATION COMPLETE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

q=3
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 
L1 = [f for f in Polq.polynomials(of_degree=1) if f.is_irreducible()
and f.leading_coefficient()==1]
L2 = [f for f in Polq.polynomials(of_degree=2) if f.is_irreducible()
and f.leading_coefficient()==1]
L3 = [f for f in Polq.polynomials(of_degree=3) if f.is_irreducible()
and f.leading_coefficient()==1]
    
f = L1[0]*L1[1]^3*L1[2]^4
f *= L2[0]*L2[1]^4*L2[2]^4
f *= L3[0]*L3[1]*L3[2]^2*L3[3]^2*L3[4]^3*L3[5]^3*L3[6]^4*L3[7]^4


# Code pour l'EXERCICE

def myFactorisation(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
    without_carre = myFsFC(f)
    for h,i in without_carre: # on obtient gi, sans carré, produit de degrés distincts
        new_h = myFEDD(h)
        for j in range(len(new_h)):
            h_fact = myCZ(new_h[j],j+1)
            if h_fact!=[1]:
                for m in h_fact:
                    retour.append((m,i))
    assert(prod([f^e for (f,e) in retour ]) == f)
    return retour


# # Affichage des resultats

print( "\n$ Question 1")
print( "La factorisation de",f,"est:\n")
for (h,i) in myFactorisation(f):
    print(f'({h})^{i}')


print("\nVersion Sage f.factor():",f.factor())


#reset()
print("""\
# ****************************************************************************
# RACINES D'UN POLYNOME
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

q=3
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 
L1 = [f for f in Polq.polynomials(of_degree=1) if f.is_irreducible()
and f.leading_coefficient()==1]
L2 = [f for f in Polq.polynomials(of_degree=2) if f.is_irreducible()
and f.leading_coefficient()==1]
L3 = [f for f in Polq.polynomials(of_degree=3) if f.is_irreducible()
and f.leading_coefficient()==1]

f = L1[0]*L1[1]^3*L1[2]^4
f *= L2[0]*L2[1]^4*L2[2]^4
f *= L3[0]*L3[1]*L3[2]^2*L3[3]^2*L3[4]^3*L3[5]^3*L3[6]^4*L3[7]^4


# Code pour l'EXERCICE


# def myRacine(f):       # version pour calculer l'entièreté des racines, y compris dans F_q^d, où d est le degré du facteur dans Fq[x]
#     Pol=f.parent()
#     x=Pol.gen()
#     p=Pol.base_ring().characteristic()
#     q=Pol.base_ring().cardinality()
#     retour = []
#     h = myFactorisation(f)
#     sorted_degree = []
#     idx = 0
#     d = 1
#     while idx<len(h):
#         if h[idx][0].degree()==d:
#             sorted_degree[d-1].append(h[idx][0])
#         else:
#             d+=1
#             sorted_degree[d-1] = [h[idx][0]]
#         idx+=1
#     for l in sorted_degree:
#         d = l[0].degree()
#         Fq_d.<alpha+str(d)>=FiniteField(q^d)
#         for p in l:
#             retour.append(fact.roots())
#     assert(f(z)==0 for z in retour)
#     return retour

def myRacine(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
    h = myFactorisation(f)
    one_degree = []
    idx = 0
    while idx<len(h):
        if h[idx][0].degree()==1: # on reste dans Fq^1[x]/(p) où deg(p)=1 car sinon obligation d'extension de corps: on sort de Fq
            one_degree.append(h[idx][0])
        idx+=1
    for l in one_degree:
        retour.append(-l[0]) # pour p = x-a, a = -p(0) car unitaires
    assert(f(z)==0 for z in retour)
    return retour

# # Affichage des resultats

print( "\n$ Question 1")
print( "Les racines de ",f,"sont",myRacine(f))



#reset()
print("""\
# ****************************************************************************
# ETUDE DE CANTOR-ZASSENHAUSS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

# Code pour l'EXERCICE
q=3
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq)
d = 5
Ld = [f for f in Polq.polynomials(of_degree=d) if f.is_irreducible() and f.leading_coefficient()==1]
f1 = random.choice(Ld)
f2 = random.choice(Ld)
while f1==f2:
    f2 = random.choice(Ld)
f = f1*f2

L1 = [g for g in Ld if (g^((q^d-1)//2)-1)%f1==0]
L1_prim = [g for g in Ld if g not in L1]
L2 = [g for g in Ld if (g^((q^d-1)//2)-1)%f2==0]
L2_prim = [g for g in Ld if g not in L2]

# # Affichage des resultats
print("\nOn a construit f = f1*f2 =",f,"avec f1 =",f1,"et f2 =",f2)

print("\nCas où u mod f1 dans L1 et u mod f2 dans L2:")
u1 = random.choice(L1)
u2 = random.choice(L2)
u = crt(u1,u2,f1,f2)
print("u est un carré dans les deux cas donc par irréductibilité et reste chinois il est à la fois carré de f1 et f2, soit u^(q^d-1/2)-1=0[f1] = 0[f2]")
print("Le gcd avec f donne f tout entier:")
print(f.gcd(exp_fast(u,(q^d-1)//2,f)-1))

print("\nCas où u mod f1 dans L1' et u mod f2 dans L2:")
u1 = random.choice(L1_prim)
u2 = random.choice(L2)
u = crt(u1,u2,f1,f2)
print("u est un carré pour f2, soit u^(q^d-1/2)-1 = 0[f2]")
print("Le gcd avec f donne f2:")
print(f.gcd(exp_fast(u,(q^d-1)//2,f)-1))

print("\n Cas où u mod f1 dans L1 et u mod f2 dans L2':")
u1 = random.choice(L1)
u2 = random.choice(L2_prim)
u = crt(u1,u2,f1,f2)
print("u est un carré pour f1, soit u^(q^d-1/2)-1 = 0[f1]")
print("Le gcd avec f donne f1:")
print(f.gcd(exp_fast(u,(q^d-1)//2,f)-1))

print("\nCas où u mod f1 dans L1' et u mod f2 dans L2':")
u1 = random.choice(L1_prim)
u2 = random.choice(L2_prim)
u = crt(u1,u2,f1,f2)
print("Aucun carré, donc u^(q^d-1)/2 = -1[fi] != 1[fi]")
print("Le gcd est trivial: 1")
print(f.gcd(exp_fast(u,(q^d-1)//2,f)-1))


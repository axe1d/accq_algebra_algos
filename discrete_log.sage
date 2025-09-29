print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP14 : LOG DISCRET ET COUPLAGES                                             #
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
# PAS DE BEBE, PAS DE GEANT
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p1 = 1823
Fp1 = FiniteField(p1)
b1 = Fp1(3)
x1 = Fp1(693)

p2 = 239
Fp2 = FiniteField(p2)
b2 = Fp2(2)
x2 = Fp2(15)


# Code pour l'EXERCICE

def Shanks(x,b):
    Fp = x.parent()
    n = Fp.order()
    s = floor(sqrt(n))+1
    T = {}
    for j in range(s):
        betaj = x*b^(-j)
        T[betaj] = j
    i = 0
    gamma = Fp(1)
    b_p = b^s
    while not (gamma in T.keys()):
        i+=1
        gamma = gamma*b_p
    j = T[gamma]
    return i*s+j


# # Affichage des resultats

print("Question 2 :", Shanks(x1,b1))
print("Question 3 :", Shanks(x2,b2))




#reset()
print("""\
# ****************************************************************************
# RHO DE POLLARD
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p= 281
Fp = FiniteField(p)
x1 = Fp(263) 
b1 = Fp(239)
x2 = Fp(165)
b2 = Fp(127)
x3 = Fp(210)
b3 = Fp(199)
x4 = Fp(2^4)
b4 = Fp(2)

# Code pour l'EXERCICE

def phi(Fp,x,b,w,alpha,beta):
    cast = lambda b: hash(b)%3
    n = Fp.order()-1
    if cast(w)==0:
        return Fp(x*w),alpha,(beta+1).mod(n)
    elif cast(w)==1:
        return Fp(w**2),(2*alpha).mod(n),(2*beta).mod(n)
    else:
        return Fp(b*w),(alpha+1).mod(n),beta

def rho(x,b):
    Fp = x.parent()
    n = Fp.order()-1
    alpha0 = Integers(n).random_element()
    beta0 = Integers(n).random_element()
    w0 = b**alpha0*x**beta0
    xx,alphax,betax = phi(Fp,x,b,w0,alpha0,beta0)
    y,alphay,betay = phi(Fp,x,b,xx,alphax,betax)
    while xx!=y:
        xx,alphax,betax = phi(Fp,x,b,xx,alphax,betax)
        y,alphay,betay = phi(Fp,x,b,*phi(Fp,x,b,y,alphay,betay))
    try:
        #print(alphax,alphay,betax,betay)
        return (alphax-alphay).mod(n)*(betay-betax)**-1
    except ZeroDivisionError:
        return rho(x,b)
        return
    #return log(x,b)


# # Affichage des resultats

print("Le log de x=",x1,"en base",b1,"vaut",rho(x1,b1),",vérification sage built-in",discrete_log(x1,b1),".")
print("Le log de x=",x2,"en base",b2,"vaut",rho(x2,b2),",vérification sage built-in",discrete_log(x2,b2),".")
print("Le log de x=",x3,"en base",b3,"vaut",rho(x3,b3),",vérification sage built-in",discrete_log(x3,b3),".")



#reset()
print("""\
# ****************************************************************************
# COUPLAGE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 61
Fp = FiniteField(p)
E = EllipticCurve(Fp,[11,0])
print("groupe de E=", E.abelian_group()) # pour verifier
S = E(24,34)
T = E(5,27)
r = 10
print("Verification de la r-torsion : r*S =", r*S, "et r*T =", r*T)


# Code pour l'EXERCICE

# def myLine(P1,P2,S):
#     E=P1.curve()
#     x1=P1[0]
#     y1=P1[1]
#     z1=P1[2]
#     x2=P2[0]
#     y2=P2[1]
#     z2=P2[2]
#     xS=S[0]
#     yS=S[1]
#     zS=S[2]

#     if P1==P2:
#         return (yS*z1-y1*zS)-(xS*z1-x1*zS)
#     else:
#         return (xS*(y1*z2-y2*z1)+yS*(x2*z1-x1*z2)+zS*(x1*y2-x2*y1))

def myLine(P1,P2,S):
    E=P1.curve()
    K = E.base_field()
    alpha = E.a4()
    x1=P1[0]; y1=P1[1]; z1=P1[2]
    x2=P2[0]; y2=P2[1]; z2=P2[2]
    xS=S[0]; yS=S[1]; zS=S[2]
    a = 1

    if z1==0 and z2==0:
        return K(1)
    elif z1==0:
        return xS-x2*zS
    elif z2==0 or (x1==x2 and y1==-y2):
        return xS-x1*zS
    elif x1==x2:
        a = (3*x1^2+alpha)/(2*y1)
    else:
        a = (y1-y2)/(x1-x2)
    return yS-y1*zS-(xS-x1*zS)*a


def myH(P1,P2,S):
    return myLine(P1,P2,S)/(myLine(P1+P2,-P1-P2,S))

def myMiller(r,S,P):
    R = S
    f = 1
    r = [int(x) for x in list(r.binary())]
    r.reverse()
    l = len(r)
    for i in range(l-2,-1,-1):
        f = myH(R,R,P)*f**2
        R = 2*R
        if r[i]==1:
            f = f*myH(R,S,P)
            R+=S
    return f

def myTatePairing(S,T,r):
    try:
        return myMiller(r,S,T)
    except ZeroDivisionError:
        Q = S.curve().random_point()
        return myTatePairing(S,T+Q,r)

def myWeilPairing(S,T,r):
    return myTatePairing(S,T,r)/myTatePairing(T,S,r)

# # Affichage des resultats

print("Calcul du couplage de Tate", myTatePairing(S,T,r), "vérification built-in Sage:", S.tate_pairing(T,r,1))
print("Calcul du couplage de Weil", myWeilPairing(S,T,r), "vérification built-in Sage:", S.weil_pairing(T,r))



#reset()
print("""\
# ****************************************************************************
# ATTAQUE M.O.V.
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 2199023255579
Fp = FiniteField(p)
E = EllipticCurve(Fp,[1,0])
P = E(1435967701832 , 123951463462)
Q = E(1129476910351 , 1383670460733)

# Code pour l'EXERCICE

j = E.j_invariant() # j-invariant a faire calculer par une fonction de SageMath
rep2 = "j=1728 et p est impair, donc p est congru à 3 modulo 4"
t = 1 # Ecrire le code pour calculer cette valeur
r = P.order()
while (p^t-1).mod(r)!=0:
    t+=1
q = p^t
Fq.<alpha> = FiniteField(q)
EE = EllipticCurve(Fq,[1,0])
PP = EE(1435967701832 , 123951463462)
QQ = EE(1129476910351 , 1383670460733)
SS = EE.random_element() # point a calculer vous-meme
rr = PP.order()
while SS.order()!=rr or PP.weil_pairing(SS,r)==1:
    SS = EE.random_element()
zeta1 = Fq(PP.weil_pairing(SS,r))
zeta2 = Fq(QQ.weil_pairing(SS,r))
lambd = log(zeta2,zeta1)


# # Affichage des resultats

print("p premier ?",p.is_prime())
print("j-invariant de E :",j)
print("p mod 4 =", mod(p,4))
print(rep2)
print("Cardinal de E(Fp) :",E.cardinality(),"=",E.cardinality().factor())
print("Ordre de P :",P.order())
print("Cardinal de E(Fq) :",EE.cardinality(),"=",EE.cardinality().factor())
print("Point S :",SS)
print("On calcule zeta1 =",zeta1,", zeta2 =",zeta2,", lambda =",lambd,".")




#reset()
print("""\
# ****************************************************************************
# CALCUL D'INDICE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 439
Fp = FiniteField(p)
g = Fp(237)
b = Fp(136)
y = 11

# Code pour l'EXERCICE

def div_successives_friable(n, P):
    f = []
    for d in P:
        cnt_d = 0
        while n.mod(d)==0:
            n = n//d
            cnt_d+=1
        f.append((d,cnt_d))
        d+=1
    if n!=1:
        return None
    return f

def prod_fri(l):
    x = 1
    for (d,cnt) in l:
        x*=d**cnt
    return x

def LogIndice(g,b,y):
    Fp = g.parent()
    p = Fp.order()
    P = prime_range(y+1)
    k = len(P)-1
    i = 1
    v = []
    a = []
    while i<=4*k:
        alpha = Integers(p-1).random_element()
        gamma = (b**alpha).mod(p)
        friy = div_successives_friable(Integer(gamma),P)
        if friy is not None and prod_fri(friy)==gamma:
            a.append(alpha)
            v.append([cnt for (_,cnt) in friy])
            i+=1
    M = matrix(v)
    e = vector(a)
    logs = M.solve_right(e)

    beta = Integers(p-1).random_element()
    beta2 = b**beta*g
    l = div_successives_friable(Integer(beta2),P)
    while l is None or prod_fri(l)!=beta2:
        beta = Integers(p-1).random_element()
        beta2 = b**beta*g
        l = div_successives_friable(Integer(beta2),P)
    f = vector([cnt for (_,cnt) in l])
    return -beta+f.dot_product(logs)

# # Affichage des resultats

print("Le log de g=",g,"en base",b,"vaut",LogIndice(g,b,y),",vérification sage built-in:",discrete_log(g,b),".")


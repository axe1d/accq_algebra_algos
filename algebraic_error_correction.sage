import random

print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP11 : CODES CORRECTEURS D'ERREURS ALGEBRIQUES                              #
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
# CODE DE HAMMING
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

C = codes.HammingCode(GF(2),3)

# Code pour l'EXERCICE

B = C.generator_matrix()
d = C.minimum_distance()

# # Affichage des resultats

print( "Le code",C,"a pour matrice")
print( B)
print( "et pour distance minimale d=",d)


#reset()
print("""\
# ****************************************************************************
# DECODAGE DE BERLEKAMP-MASSEY
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

from sage.matrix.berlekamp_massey import berlekamp_massey

q=13
Fq = FiniteField(q)
k=5
alpha=Fq(6)
PolFq.<x> = PolynomialRing(Fq)
r = vector(Fq,[5,2,9,9,1,11,7,5,7,4,3,1])
#r = vector(Fq,[0,0,5,4,0,4,2,0,7,7])
n = q-1

# Code pour l'EXERCICE

G = matrix(Fq,[[z^j for z in [alpha^i for i in range(n)]] for j in range(k)])
H = matrix(Fq,[[z^j for z in [alpha^i for i in range(n)]] for j in range(1,n-k+1)])
s = H*r

if len(list(s))%2==1:
    s = vector(Fq,[z for i,z in enumerate(list(s)) if i!=len(list(s))-1])
sigma = berlekamp_massey(list(s))
#sigma,_ = myBerlekampMassey(list(s))
M = [z.log(alpha) for z in sigma.roots(multiplicities=False)]
r_new = vector(Fq,[z for (i,z) in enumerate(list(r)) if i not in M])
G_new = matrix(Fq,[[z^j for z in [alpha^i for i in range(n) if i not in M]] for j in range(k)])
m = G_new.solve_left(r_new)

# # Affichage des resultats

print( "La matrice génératrice du code est")
print( G)
print( "La matrice de controle du code est")
print( H)
print( "Le syndrome est")
print( s)
print( "Le polynome localisateur d'erreurs est")
print( sigma)
print( "La position des erreurs est")
print( M)
print( "Le message envoye le plus probable est")
print( m)



#reset()
print("""\
# ****************************************************************************
# DECODAGE EN LISTE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

q=23
#q=17
Fq = FiniteField(q)
k = 3

alpha=Fq(14)
#alpha=Fq(6)
MPol.<x,y>=PolynomialRing(Fq,2)

r=vector(Fq, [12,18,15,22,17,5,14,21,17,4,13,8,4,10, 15,11,22,12,13,9,14,12])
#r = vector(Fq,[15 ,16 ,7 ,9 ,16 ,6 ,15 ,0 ,7 ,7 ,4 ,8 ,9 ,5 ,3 ,6])
n = q-1

G = matrix(Fq,[[z^j for z in [alpha^i for i in range(n)]] for j in range(k)])

Rxy = PolynomialRing(Fq,['x','y'])
x,y = Rxy.gens()

# Code pour l'EXERCICE

b1=floor((n-k)/2)
b2=floor(n-sqrt(2*(k-1)*n))

cand_mono = []
for i in range(b2):
    for j in range(b2):
        if i+(k-1)*j<b2:
            cand_mono.append(x^i*y^j)

lignes = [(alpha^i,list(r)[i]) for i in range(n)]
mat_Q = matrix(Fq,[[z(x,y) for z in cand_mono] for (x,y) in lignes])
K = mat_Q.right_kernel().basis()
list_poly_gen = []
for p in K:
    pol = MPol(0)
    for (i,z) in enumerate(list(p)):
        pol+=z*cand_mono[i]
    list_poly_gen.append(pol)


fact_candidat = []
for p in list_poly_gen:
    irred = True
    for (c,_) in list(p.factor()):
        if c.degree(y)==1:
            irred = False
    if not irred:
        fact_candidat.append(p)


print(fact_candidat[0])
list_fi = []
for (p,_) in list(factor(fact_candidat[0])):
    if p.mod(y)!=p:
        list_fi.append(y-p)

list_fi=list_fi[:2]

L = [vector(Fq,[f(alpha^i,0) for i in range(n)]) for f in list_fi]
M = [G.solve_left(vector(Fq,list(c))) for c in L]
M = [(k,[chr(ord('A')+int(j)) for j in list(k)]) for k in M]


def hamming(u, v):
    return sum(1 for i in range(len(u)) if u[i] != v[i])
d1=hamming(L[0],r)
d2=hamming(L[1],r)

c = L[0] if d1<d2 else L[1]
m = M[0][0] if d1<d2 else M[1][1]

# # Affichage des resultats

print( "Le decodage unique fonctionne jusqu'à distance",b1)
print( "Le decodage de Sudan fonctionne jusqu'à distance",b2)
print( "La liste des mots de code les plus proches de r est",L)
print( "Elle correspond aux messages")
print( M)
print( "Le mot de code le plus probable est",c)
print( "soit le message",m)


#reset()
print("""\
# ****************************************************************************
# CODE ALGEBRIQUE
# ****************************************************************************
""")


# ESSAYEZ DE LE TRAITER
# CET EXERCICE NE SERA PAS EVALUE
# IL VOUS MANQUE LE THEOREME DE RIEMANN-ROCH POUR TOUT COMPRENDRE


# q  = 23
# Fq = GF(q)
# # Example: y^2 = x^3 + a*x + b over Fq
# a, b = Fq(1), Fq(1)
# E  = EllipticCurve(Fq, [a,b])
# O  = E(0)

# raw_pts = [P for P in E.points() if P != O]
# Pts     = [(P[0], P[1]) for P in raw_pts]   # now a list of 2-tuples (x,y)
# n   = len(Pts)
# print(f"Curve has n={n} rational affine points.")

# m = 5

# # 4. Build the L(mO)-basis as monomials x^i y^j with 2i+3j ≤ m
# Rxy = PolynomialRing(Fq,'x,y')
# x,y = Rxy.gens()
# monomials_LG = [ x^i * y^j 
#                  for i in range(m+1) 
#                  for j in range(m+1)
#                  if 2*i + 3*j <= m ]

# k = len(monomials_LG)
# print("dim L(G) =", k)

# # 5. Generator matrix G: k×n
# G_ec = matrix(Fq, k, n,
#     lambda r,c: monomials_LG[r]( *Pts[c] ) )


# def encode_ec(m):
#     return vector(Fq,m) * G_ec

# # 6. List‐decoding radius b2 (for Sudan)
# b2 = floor(n - sqrt(k*n))      # e.g. Sudan bound

# # 7. Candidate monomials for Q(x,y):
# #    choose all x^i y^j with i + (k-1)*j < b2
# cand_Q = [ x^i * y^j
#            for i in range(b2)
#            for j in range(b2)
#            if i + (k-1)*j < b2 ]

# # 8. Build the evaluation matrix of Q at (P_i, r_i)
# #    Here r is your received vector of length n.
# xPts = [ P[0] for P in E.points() if P != O ]
# make_Q_matrix = lambda r: matrix(Fq, n, len(cand_Q),
#     lambda i,j: cand_Q[j]( xPts[i], r[i] ) )

# # def make_Q_matrix(r):
# #     return matrix(Fq, n, len(cand_Q),
# #        lambda i,j: cand_Q[j](*Pts[i], r[i]))

# # 9. Kernel → list of Q’s
# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# MPol = PolynomialRing(Fq, ['x','y'])
# def list_candidates(r):
#     Qmat = make_Q_matrix(r)
#     K    = Qmat.right_kernel().basis()
#     polys = []
#     for row in K:
#         polys.append(sum(row[j]*cand_Q[j] for j in range(len(cand_Q))))
#     return polys

# # 10. Factor each Q, pull out any factor y - f(x):
# def extract_f(polys):
#     Fs = []
#     for Q in polys:
#         for F,mul in Q.factor():
#             if F.degree(y)==1:
#                 # F = a(x)*y + b(x)
#                 a = F.coefficient({x:0,y:1})
#                 b = F.subs(y=0)
#                 f = -b/a
#                 Fs.append(f)
#     return Fs

# # 11. From each f(x) ∈ L(G) you get a candidate codeword:
# def list_decode(r):
#     polys = list_candidates(r)
#     Fs    = extract_f(polys)
#     Cs    = [ vector(Fq, [fi(*Pts[i]) for i in range(n)]) for fi in Fs ]
#     return Cs

# def hamming(u,v):
#     return sum(1 for i in range(n) if u[i]!=v[i])

# def decode_ec(r):
#     Cs = list_decode(r)
#     return min(Cs, key=lambda c: hamming(c,r))
# # Example: send a random message
# msg = [ Fq.random_element() for _ in range(k) ]
# cw  = encode_ec(msg)

# # introduce ≤ ⌊(b2-1)/2⌋ errors
# r   = vector(Fq,cw)
# for idx in random.sample(range(n), 3):
#     r[idx] += Fq.random_element()

# # decode
# ĉ = decode_ec(r)
# print("Recovered codeword:", ĉ)
# print("Original was      :", cw)




x,y = polygens(GF(2),'x, y')
Pol.<x,y>=PolynomialRing(GF(2), 'x,y')
F2=GF(2)
F4.<a>=GF(4,'a')
F8.<b>=GF(8,'b')
F16.<c>=GF(16,'c')
p = y^2+y+x^3+x+1
ptsaffines4=[ (u,v) for u in F4 for v in F4 if p(u,v)==0]
ptsaffines8=[ (u,v) for u in F8 for v in F8 if p(u,v)==0]
8.
ptsaffines16=[ (u,v) for u in F16 for v in F16 if p(u,v)==0]
print(len(ptsaffines4), len(ptsaffines8), len(ptsaffines16))
L_Base=[Pol(1)]
for i in range(1,5):
    L_Base.append(x^i)
    L_Base.append(x^(i-1)*y)
L_Code = [[ f(u) for u in ptsaffines8] for f in L_Base[0:7]]
Code = LinearCode(Matrix(L_Code))
print(Code.length())
print(Code.dimension())
print(Code.minimum_distance())



#reset()
print("""\
# ****************************************************************************
# BORNES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

#q=49

# Code pour l'EXERCICE
entropy_q = lambda q: lambda x: 0 if x in [0,1] else (x*log(q-1,q)-x*log(x,q)-(1-x)*log(1-x,q))
theta = lambda q: (q-1)/q

gv = lambda q: lambda x: 1 if (x==0 or x>1+1/q) else 1-entropy_q(q)(x)
plotkin = lambda q: lambda x: 0 if x>=theta(q) else 1-x/theta(q)
hamming = lambda q: lambda x: 0 if x/2>theta(q) else 1-entropy_q(q)(x/2)
mrrw = lambda q: lambda x: 0 if x>theta(q) else entropy_q(q)((1/q)*(q-1-(q-2)*x-2*sqrt((q-1)*x*(1-x))))
better_than_gv = lambda q: lambda x: 1-x-1/(sqrt(q)-1)

# # Affichage des resultats

def plot_bounds(q):
    plot([gv(q),plotkin(q),hamming(q),mrrw(q),better_than_gv(q)],xmin=0,xmax=theta(q),ymin=0,ymax=1, legend_label=["GV","Plotkin","Hamming","MRRW","TVF>GV"],title="Differentes bornes sur les codes").save('code-bounds-q'+str(q)+'.png')

#plot_bounds(49)
#plot_bounds(200)
#plot_bounds(1000)

print("""\
# ****************************************************************************
# BORNES BATTANT GILBERT-VARSHAMOV
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

# Code pour l'EXERCICE

# codesalg = .7

# # # Affichage des resultats

# plot([gv,plotkin,hamming,mrrw,codesalg],xmin=0,xmax=1,ymin=0,ymax=1, legend_label=["GV","Plotkin","Hamming","MRRW","Codes Algebriques"],title="Differentes bornes sur les codes")

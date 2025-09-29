print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP12 : CRYPTANALYSE ET CRYPTOGRAPHIE A BASE DE RESEAUX                      #
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
# SAC A DOS (NE PAS TRAITER LA QUESTION 2)
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

b = [356,278,417,27,132,464,521]
s = 1287

# Code pour l'EXERCICE

B = matrix.identity(len(b))
B = B.insert_row(len(b),list(b))
B = B.transpose()
B = B.insert_row(len(b),[0]*len(b)+[-s])

x = B.LLL()[0]
x = vector(list(x)[:-1])

# # Affichage des resultats

print("Le message est")
print(x)

#reset()
print("""\
# ****************************************************************************
# ATTAQUE DE WIENER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

Pol.<x> = PolynomialRing(ZZ)

N1 = 65946239999
e1 = 22022476093
#N1 = 24585612803
#e1 = 10878450977
N2 = 65946239999
e2 = 10865199773

# Code pour l'EXERCICE
def find_p_q(N1,e1):
    X1 = round(N1**(1/4))
    Y1 = round(N1**(3/4))
    b1 = vector(ZZ,[e1*X1,Y1])
    b2 = vector(ZZ,[N1*X1,0])
    M = matrix(ZZ,[b1,b2])
    v = M.LLL().row(0)
    coeffs = M.solve_left(v)
    k = abs(coeffs[1])
    s = abs(coeffs[0])
    sigma = (k*N1-e1*s+1)//k+1
    f = x^2-sigma*x+N1
    roots = f.roots(multiplicities=False)
    if roots!=[]:
        p = roots[1]
        q = roots[0]
        print("Clé cassée")
        return p,q
    else:
        print("Aucune solution trouvée: bonne clé!")

print("Première clé",find_p_q(N1,e1))
print("\n")
print("Deuxième clé",find_p_q(N2,e2))

ans  ="1"

# # Affichage des resultats

print("Il vaut mieux utiliser la cle ",ans)

#reset()
print("""\
# ****************************************************************************
# METHODE DE COPPERSMITH
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

Pol.<x> = PolynomialRing(ZZ)

# Code pour l'EXERCICE

def Coppersmith(f,N,B=None):
    R = f.parent()
    x = R.gen()
    d = f.degree()
    m = ceil(log(N)/d)
    if B is None:
        B = ceil(1/(2*exp(1))*N**(1/d))

    M = []
    l_g_ij = []
    for i in range(m+1):
        for j in range(d):
            gij = expand(x^j*N^i*f^(m-i))
            gij = list(gij(x=B*x))
            M.append(gij+[0]*((m+1)*d-len(gij)))
    M = matrix(ZZ,M)
    v = M.LLL().row(0)

    h = R(list(v))
    h = R([k//B^i for (i,k) in enumerate(list(h))])

    R = h.roots(multiplicities=False)
    return R

# # Affichage des resultats

p=(x+1)*(x-2)*(x-3)*(x-29)
print(p,"\nPetites racines avec Coppersmith:",Coppersmith(p,10000))
p = x^2-27*x-33
print(p,"\nPetites racines avec Coppersmith:",Coppersmith(p,143))


#reset()
print("""\
# ****************************************************************************
# MESSAGES STEREOTYPES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


bin=BinaryStrings()
N = 42564360034887861127
Pol.<x> = PolynomialRing(ZZ)
PolmodN.<y> = PolynomialRing(Integers(N))
e = 3
text = "08/06:"
m = ZZ(str(bin.encoding(text)), base=2)*2^16

c = 12843085802751039909
f = expand((m+x)**e-c)

#r = Coppersmith(f,N,B=2^16)[0]
r = 15411
print("Racine trouvée par Coppersmith:",r)

# Code pour l'EXERCICE

#mm=ZZ(72*256^3 + 69*256^2+ 76*256 + 80).str(2) # A CHANGER DANS LES ()
mm=r.str(2)
mm = '0'*(8*ceil(len(mm)/8)-len(mm))+mm
mm = bin(mm).decoding()

# # Affichage des resultats

print("Ce jour la, le message est")
print(mm)


#reset()
print("""\
# ****************************************************************************
# ALGORITHME DE BABAI (NE TRAITER QUE LA QUESTION 1)
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

def Babai(B,t):
    b = t
    G,_ = B.gram_schmidt()
    for j in range(B.nrows()-1,-1,-1):
        uj = b.dot_product(G.row(j))/(G.row(j).norm())^2
        b = b-round(uj)*B.row(j)
    return t-b

# # Affichage des resultats


#reset()
print("""\
# ****************************************************************************
# CRYPTOSYSTEME GGH
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

EE= matrix(ZZ,40, {(32, 32): 1, (22, 13): -2, (23, 8): -3,
       (15, 31): -1, (22, 37): -1,
       (13, 5): -4, (38, 20): 2, (4, 12): 3, (19, 22): -2, (15, 5): 2, (11,
       32): -1, (11, 10): 3, (1, 11): -4, (12, 33): 1, (0, 15): 1, (33, 17): 1,
       (7, 19): -1, (11, 1): -2, (7, 27): 3, (19, 32): -4, (22, 10): 2, (31,
       39): -4, (34, 9): 2, (36, 17): 2, (18, 17): 1, (14, 6): -2, (23, 14): 3,
       (23, 34): 2, (12, 11): -3, (0, 21): -3, (27, 22): -2, (4, 29): -3, (23,
       5): 1, (4, 6): -2, (24, 7): 2, (5, 38): -2, (33, 13): -1, (9, 35): 3,
       (18, 36): 1, (22, 5): 1, (24, 25): 3, (34, 31): 2, (6, 34): -3, (23,
       33): -4, (20, 37): -1, (38, 12): 2, (33, 0): -1, (4, 32): 3})
AA=10*identity_matrix(40)+EE
HH,UU=AA.hermite_form(transformation=True)
cc = vector([-2, 0, 2, 0, 0, 1, -1, -1, -3, 0, 0, 2, -1, 13, 7, 2, 0, 2, 27, 2, 1,
       17, -2, 899, 50, 15, 11, 1098, 7, 2, -1, 10, -1, 2, 156, 15, 42, 8,
       525748584, 37])


# Code pour l'EXERCICE
m2 = Babai(AA,cc)
m = m2*HH.inverse()

def StringToAscii(text):
    return [ZZ(ord(s)) for s in list(text)]
def AsciiToString(listascii):
    return "".join([str(unichr(a)) for a in listascii])
bin = BinaryStrings()
def StringToInts(text):
    return [ZZ(str(t),base=2) for t in list(bin.encoding(text))]
def IntsToString(listints):
    return bin("".join([str(t) for t in listints])).decoding()

print(IntsToString(list(m)))

# # Affichage des resultats


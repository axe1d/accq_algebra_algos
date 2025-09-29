print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP2 : ALGEBRE LINEAIRE SUR UN ANNEAU PRINCIPAL                              #
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
# MISE SOUS FORME NORMALE D'HERMITE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
        [-2,  3,  3,  1],
        [ 2, -1,  1, -3],
        [-4,  0, -1, -4]])

A1 = random_matrix(ZZ, 7, 8, algorithm='echelonizable', rank=3)

U = identity_matrix(4)

# Code pour l'EXERCICE

def cherche_pivot_non_nul(A,U,i,j):
    """
    Echange la colonne j avec une colonne à gauche pour que A[i,j] soit non nul.
    """
    if i<0 or i>=A.nrows() or j<0 or j>=A.ncols():
        return False, A, U
    r = j
    while r >= 0 and A[i, r] == 0:
        r -= 1
    if r < 0 or A[i, r] == 0:
        return False,A,U
    if r!=j:
        cj = A[:,j]
        A[:,j] = A[:,r]
        A[:,r] = cj
        U = U*elementary_matrix(U.ncols(),row1=r,row2=j)
    return True,A,U

def normalise_pivot(A,U,i,j):
    """
    Multiplie la colonne j si besoin pour que A[i,j] soit positif.
    """
    assert(A[i,j]!=0)
    if A[i,j]<0:
        A[:,j] *= -1
        U = U*elementary_matrix(U.ncols(),row1=j,scale=-1)
    return A,U

def annule_a_gauche(A,U,i,j):
    """
    Annule les coefficients à gauche de A[i,j]
    """
    assert(A[i,j]!=0)
    for r in range(j-1,-1,-1):
        if A[i,r]!=0:
            # on échange les colonnes pour appliquer Bezout
            cj = A[:,j]
            A[:,j] = A[:,r]
            A[:,r] = cj
            U = U*elementary_matrix(U.ncols(),row1=r,row2=j)
            d,u,v = xgcd(A[i,r],A[i,j])
            if d == 0:
                continue
            # stockage en mémoire valeurs
            ci = A[:,r]
            ciu = U[:,r]
            cj = A[:,j]
            cju = U[:,j]
            a = A[i,r]
            b = A[i,j]

            A[:,r] = u*ci # Ci <- uCi
            U[:, r] = u * ciu
            A[:,r] = A[:,r] + v*cj # Ci <- Ci+vCj
            U[:,r] = U[:,r] + v*cju
            A[:,j] = (a/d)*cj -(b/d)*ci # Cj <- (a/d)Cj -(b/d)Ci
            U[:,j] = (a/d)*cju -(b/d)*ciu
            A = A*elementary_matrix(A.ncols(),row1=r,row2=j) # on échange à nouveau
            U = U*elementary_matrix(U.ncols(),row1=r,row2=j)
    return A,U


def reduit_a_droite(A,U,i,j):
    """
    Réduit les coefficients à gauche de A[i,j] modulo A[i,j]
    """
    assert(A[i,j]!=0)
    for r in range(j+1,A.ncols()):
        if A[i, r] != 0:
            d = A[i, r]//A[i,j]
            A[:, r] -= d*A[:,j]
            U[:, r] -= d*U[:,j]
    return A,U


def MyHNF(A):
    """
    Forme normale d'Hermite selon votre code
    """
    m = A.nrows()
    n = A.ncols()
    H = copy(A)
    U = identity_matrix(n)
    # ECRIVEZ VOTRE CODE ICI, VOUS POUVEZ REPRENDRE LES FONCTIONS PRECEDENTES
    # COMME SOUS-FONCTION
    for j in range(n-1,-1,-1):
        i = m-1
        has_pivot, H, U = cherche_pivot_non_nul(H,U,i,j)
        while not has_pivot and i>=0:
            i -= 1
            if i>=0:
                has_pivot, H, U = cherche_pivot_non_nul(H,U,i,j)
        if not has_pivot:
            continue
        H, U = annule_a_gauche(H,U,i,j)
        H, U = normalise_pivot(H,U,i,j)
        H, U = reduit_a_droite(H,U,i,j)
    assert(H-A*U==0)
    return H,U


def SageHNF(A):
    """
    Forme normale d'Hermite de SAGE avec la convention francaise :
    Les vecteurs sont ecrits en colonne.
    """
    m = A.nrows()
    n = A.ncols()
    Mm = identity_matrix(ZZ,m)[::-1,:]
    Mn = identity_matrix(ZZ,n)[::-1,:]
    AA = (Mm*A).transpose()
    HH, UU = AA.hermite_form(transformation=True)
    H = (HH*Mm).transpose()*Mn
    U = UU.transpose()*Mn
    assert(H-A*U==0)
    return H,U


H,U = MyHNF(A)
HH,UU = SageHNF(A)

A2 = [random_matrix(ZZ,ZZ.random_element(20),ZZ.random_element(20)) for i in range(100)]
test = []
for mat in A2:
    H1,U1 = MyHNF(mat)
    HH1,UU1 = SageHNF(mat)
    test.append(H1==HH1)
test = all(test)

# # Affichage des resultats

print("\n$ Question 4")
print("La matrice A = ")
print(A)
print("a pour forme normale d'Hermite H=")
print(H)
print("et matrice de transformation U=")
print(U)
print("\n$ Question 5")
print("D'apres SageMath, la matrice A a pour forme normale d'Hermite H=")
print(HH)
print("et matrice de transformation U=")
print(UU)
print("\n$ Question 6")
print("Les deux fonctions coincident-elles sur une centaine d'exemples ?")
print(test)


# Exercice 59
print("\n Exercice 59")
A = matrix(ZZ,[
        [-2,  3,  3,  1],
        [ 2, -1,  1, -3],
        [-4,  0, -1, -4]])
H,U = MyHNF(A)
print(H,U)

print("""\
# ****************************************************************************
# RESOLUTION DE SYSTEME LINEAIRE HOMOGENE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X = matrix(ZZ,[
      [ -2,  3,  3,  1],
      [  2, -1,  1, -3],
      [ -4,  0, -1, -4]])

# Code pour l'EXERCICE
H,U = MyHNF(X)
L =[vector(U[:,i]) for i in range(4) if vector(H[:,i]) == vector([0,0,0])] # liste des vecteurs d'une base

# # Affichage des resultats

print("Le systeme a pour racine le module engendre par")
print(L)

#reset()
print("""\
# ****************************************************************************
# MISE SOUS FORME NORMALE DE SMITH
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X1 = matrix(ZZ,2,3,[
        [4, 7, 2],
        [2, 4, 6]])

X2 = matrix(ZZ,3,3,[
        [-397, 423, 352],
        [   2,  -3,   1],
        [-146, 156, 128],
])

PolQ.<xQ> = PolynomialRing(QQ)
AQ = matrix(PolQ,3,[
            [xQ + 1,  2,     -6],
            [     1, xQ,     -3],
            [     1,  1, xQ - 4]])
Pol2.<x2> = PolynomialRing(FiniteField(2))
AF2 = matrix(Pol2,3,[
            [x2 + 1,  2,     -6],
            [     1, x2,     -3],
            [     1,  1, x2 - 4]])
Pol3.<x3> = PolynomialRing(FiniteField(3))
AF3 = matrix(Pol3,3,[
            [x3 + 1,  2,     -6],
            [     1, x3,     -3],
            [     1,  1, x3 - 4]])
Pol5.<x5> = PolynomialRing(FiniteField(5))
AF5 = matrix(Pol5,3,[
            [x5 + 1,  2,     -6],
            [     1, x5,     -3],
            [     1,  1, x5 - 4]])

# Code pour l'EXERCICE
# On identifie un polynôme à sa représentation matricielle dans une base du corps considéré
def cherche_pivot_non_nul_pol(E,A,U,i,j):
    """
    Echange la colonne j avec une colonne à gauche pour que A[i,j] soit non nul.
    """
    if i<0 or i>=A.nrows() or j<0 or j>=A.ncols():
        return False, A, U
    r = j
    while r >= 0 and A[i, r] == 0:
        r -= 1
    if r < 0 or A[i, r] == 0:
        return False,A,U
    if r!=j:
        cj = A[:,j]
        A[:,j] = A[:,r]
        A[:,r] = cj
        U = U*elementary_matrix(E,U.ncols(),row1=r,row2=j)
    return True,A,U

def normalise_pivot_pol(E,A,U,i,j):
    """
    Multiplie la colonne j si besoin pour que A[i,j] soit positif.
    """
    assert(A[i,j]!=0)
    if (A[i,j].coefficients())[-1]<0:
        A[:,j] *= -1
        U = U*elementary_matrix(E,U.ncols(),row1=j,scale=-1)
    return A,U

def annule_a_gauche_pol(E,A,U,i,j):
    """
    Annule les coefficients à gauche de A[i,j]
    """
    assert(A[i,j]!=0)
    for r in range(j-1,-1,-1):
        if A[i,r]!=0:
            # on échange les colonnes pour appliquer Bezout
            cj = A[:,j]
            A[:,j] = A[:,r]
            A[:,r] = cj
            U = U*elementary_matrix(E,U.ncols(),row1=r,row2=j)
            d,u,v = xgcd(A[i,r],A[i,j])
            if d == 0:
                continue
            # stockage en mémoire valeurs
            ci = A[:,r]
            ciu = U[:,r]
            cj = A[:,j]
            cju = U[:,j]
            a = A[i,r]
            b = A[i,j]

            A[:,r] = u*ci # Ci <- uCi
            U[:, r] = u * ciu
            A[:,r] = A[:,r] + v*cj # Ci <- Ci+vCj
            U[:,r] = U[:,r] + v*cju
            A[:,j] = (a/d)*cj -(b/d)*ci # Cj <- (a/d)Cj -(b/d)Ci
            U[:,j] = (a/d)*cju -(b/d)*ciu
            A = A*elementary_matrix(E,A.ncols(),row1=r,row2=j) # on échange à nouveau
            U = U*elementary_matrix(E,U.ncols(),row1=r,row2=j)
    return A,U


def reduit_a_droite_pol(A,U,i,j):
    """
    Réduit les coefficients à gauche de A[i,j] modulo A[i,j]
    """
    assert(A[i,j]!=0)
    for r in range(j+1,A.ncols()):
        if A[i, r] != 0:
            d = A[i, r]//A[i,j]
            A[:, r] -= d*A[:,j]
            U[:, r] -= d*U[:,j]
    return A,U

def myHNFforPolynomials(A,E):
    """
    Forme normale d'Hermite selon votre code
    """
    m = A.nrows()
    n = A.ncols()
    H = copy(A)
    U = identity_matrix(A.base_ring(),n)
    # ECRIVEZ VOTRE CODE ICI, VOUS POUVEZ REPRENDRE LES FONCTIONS PRECEDENTES
    # COMME SOUS-FONCTION
    for j in range(n-1,-1,-1):
        i = m-1
        has_pivot, H, U = cherche_pivot_non_nul_pol(E,H,U,i,j)
        while not has_pivot and i>=0:
            i -= 1
            if i>=0:
                has_pivot, H, U = cherche_pivot_non_nul_pol(E,H,U,i,j)
        if not has_pivot:
            continue
        H, U = annule_a_gauche_pol(E,H,U,i,j)
        H, U = normalise_pivot_pol(E,H,U,i,j)
        H, U = reduit_a_droite_pol(H,U,i,j)
    assert(H-A*U==0)
    return H,U



######################################################################################
# Forme normale de Smith
def cherche_pivot_non_nul(A,L,U,k):
    """
    Echange la ligne et colonne i et j pour que A[k,k] soit non nul.
    """
    m = A.nrows()
    n = A.ncols()
    r1 = k
    r2 = k
    while r1<m and r2<n and A[r1,r2]==0:
        if r2==n-1:
            r2 = k
            r1 += 1
        else:
            r2 += 1
    if r1>=m or A[r1, r2] == 0:
        return False,A,L,U
    if r1!=k or r2!=k:
        A[k,:],A[r1,:] = A[r1,:],A[k,:]
        A[:,k],A[:,r2] = A[:,r2],A[:,k]
        L[k,:],L[r1,:] = L[r1,:],L[k,:]
        U[:,k],U[:,r2] = U[:,r2],U[:,k]
    return True,A,L,U

def bezout_smith_i(A,L,U,k,i):
    d,s,t = xgcd(A[k,k],A[i,k])
    if d == 0 or A[i,k] == 0 or A[k,k] == 0:
        return A,L,U
    if A[i,k]%A[k,k]==0:
        d = A[k,k]
        s = 1
        t = 0
    u = (-A[i,k])/d
    v = A[k,k]/d
    li = A[i,:]
    lil = L[i,:]
    lk = A[k,:]
    lkl = L[k,:]
    A[k,:] = s*lk + t*li
    L[k,:] = s*lkl + t*lil
    A[i,:] = u*lk + v*li
    L[i,:] = u*lkl + v*lil
    return A,L,U


def bezout_smith_j(A,L,U,k,j):
    d,s,t = xgcd(A[k,k],A[k,j])
    if d == 0 or A[k,j] == 0 or A[k,k] == 0:
        return A,L,U
    if A[k,j]%A[k,k]==0:
        d = A[k,k]
        s = 1
        t = 0
    u = (-A[k,j])/d
    v = A[k,k]/d
    cj = A[:,j]
    ck = A[:,k]
    cju = U[:,j]
    cku = U[:,k]
    A[:,k] = s*ck + t*cj
    A[:,j] = u*ck + v*cj
    U[:,k] = s*cku + t*cju
    U[:,j] = u*cku + v*cju
    return A,L,U


def cond_null(A,k):
    m = A.nrows()
    n = A.ncols()
    for i in range(k+1,m):
        if A[i,k]!=0:
            return False
    for j in range(k+1,n):
        if A[k,j]!=0:
            return False
    return True

def try_not_div(A,L,U,k):
    m = A.nrows()
    n = A.ncols()
    for i in range(k+1,m):
        for j in range(k+1,n):
            if A[i,j]%A[k,k]!=0:
                A[:,k] += A[:,j]
                U[:,k] += U[:,j]
    return A,L,U

def test_zero_i(A,k):
    m = A.nrows()
    for i in range(k+1,m):
        if A[i,k]!=0:
            return False
    return True

def test_zero_j(A,k):
    n = A.ncols()
    for j in range(k+1,n):
        if A[k,j]!=0:
            return False
    return True

def normalise_pivot(A,L,U,k):
    if A[k,k]<0:
        A[:,k] *= -1
        U[:,k] *= -1
    return A,L,U


def MySNF(A):
    """
    Forme normale de Smith selon votre code
    """
    m = A.nrows()
    n = A.ncols()
    H = copy(A)
    L = identity_matrix(A.base_ring(),m)
    U = identity_matrix(A.base_ring(),n)
    # ECRIVEZ VOTRE CODE ICI
    for k in range(min(m,n)):
        H,L,U = normalise_pivot(H,L,U,k)
        while not cond_null(H,k):
            if H[k,k]==0:
                _,H,L,U = cherche_pivot_non_nul(H,L,U,k)
            H,L,U = normalise_pivot(H,L,U,k)

            for i in range(k+1,m):
                H,L,U = bezout_smith_i(H,L,U,k,i)

            if cond_null(H,k): # pour éviter une boucle infinie pour rang 1 car opération sur une colonne annule l'opération faite par la ligne
                break

            for j in range(k+1,n):
                H,L,U = bezout_smith_j(H,L,U,k,j)
            if H[k,k]!=0:
                H,L,U = try_not_div(H,L,U,k)
    assert(H-L*A*U==0)
    return H,L,U

H1, L1, U1 = MySNF(X1)
H2, L2, U2 = MySNF(X2)

HQ, _, _ = MySNF(AQ)
HF2, _, _ = MySNF(AF2)
HF3, _, _ = MySNF(AF3)
HF5, _, _ = MySNF(AF5)

A2 = [random_matrix(ZZ,ZZ.random_element(20),ZZ.random_element(20)) for i in range(100)]
test = []
for mat in A2:
    Htest,Ltest,Utest = MySNF(mat)
    HHtest = mat.smith_form()
    test.append(Htest==HHtest[0])
test = all(test)

# # Affichage des resultats

print("\n$ Question 4")
print("La matrice X1 = ")
print(X1)
print("a pour forme normale de Smith H1=")
print(H1)
print("et matrice de transformation L1=")
print(L1)
print("et matrice de transformation U1=")
print(U1)
print("La matrice X2 = ")
print(X2)
print("a pour forme normale de Smith H2=")
print(H2)
print("et matrice de transformation L2=")
print(L2)
print("et matrice de transformation U2=")
print(U2)

print("\n$ Question 5")
print("La forme normale de Smith sur Q est ")
print(HQ)
print("La forme normale de Smith sur F2 est ")
print(HF2)
print("La forme normale de Smith sur F3 est ")
print(HF3)
print("La forme normale de Smith sur F5 est ")
print(HF5)

print("\n$ Question 6")
print("Votre fonction coincide avec celle de Sage ?")
print(test)


# Question 5
print("\n Question 5")
F2 = FiniteField(2)
PolF2.<x> = PolynomialRing(F2)
A = matrix(PolF2,3,3,[
    [x+1,2,-6],
    [1,x,-3],
    [1,1,x-4]])
H,L,U = MySNF(A)
print(H)
F3 = FiniteField(3)
PolF3.<x> = PolynomialRing(F3)
A = matrix(PolF3,3,3,[
    [x+1,2,-6],
    [1,x,-3],
    [1,1,x-4]])
H,L,U = MySNF(A)
print(H)
F5 = FiniteField(5)
PolF5.<x> = PolynomialRing(F5)
A = matrix(PolF5,3,3,[
    [x+1,2,-6],
    [1,x,-3],
    [1,1,x-4]])
H,L,U = MySNF(A)
print(H)


#reset()


#reset()
print("""\
# ****************************************************************************
# IMAGE D'UNE MATRICE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A  = matrix(ZZ, [
           [ 15,  8, -9, 23,  -9],
           [ 22, 22,  7, -8,  20],
           [ 21, 18, -1, -7,  -3],
           [  3, -1,  0, 12, -16]])


# Code pour l'EXERCICE
Z4 = ZZ^4
M = Z4.submodule(A.transpose())
test = M==Z4

# # Affichage des resultats

print("L'image de")
print(A)
print("est-elle egale a ZZ^4 ?")
print(test)



#reset()
print("""\
# ****************************************************************************
# RESOLUTION DE SYSTEME LINEAIRE NON-HOMOGENE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X1 = matrix(ZZ,[
           [ -6,  12,  6],
           [ 12, -16, -2],
           [ 12, -16, -2]])

b1 = vector(ZZ,[ -6, 4, 4])

PolF5.<x> = PolynomialRing(GF(5))

X2 = matrix(PolF5,[
           [ x + 1, 2,     4],
           [     1, x,     2],
           [     1, 1, x + 1]])

b2 = vector(PolF5,[ 3*x+2, 0, -1])

F5 = FiniteField(5)

X3 = matrix(ZZ,[
           [ 2,  -3,  6, 0],
           [ 4, 4, 53, -5]])

b3 = vector(F5,[-4,2])

# Code pour l'EXERCICE

H,L,U = MySNF(X1)
b11 = L*b1
non_nul = [j for j in range(min(X1.nrows(),X1.ncols())) if H[j,j]!=0]
nul = [j for j in range(min(X1.nrows(),X1.ncols())) if H[j,j]==0]
sol_part = []
for k in non_nul:
    sol_part.append((b11[k]//H[k,k])*U[:,k])
z1 = sum(sol_part)
homogen = [vector(U[:,k]) for k in nul]
mod = span(homogen,ZZ) # module de la solution homogène
H1 = [z1,mod]

H,L,U = MySNF(X2)
print(H,L,U)
b22 = L*b2
print(b2,b22)
non_nul = [j for j in range(min(X2.nrows(),X2.ncols())) if H[j,j]!=0]
nul = [j for j in range(min(X2.nrows(),X2.ncols())) if H[j,j]==0]
sol_part = []
for k in non_nul:
    q= b22[k]/(H[k,k])
    sol_part.append(q*U[:,k])
z2 = sum(sol_part)
homogen = [vector(U[:,k]) for k in nul]
mod2 = span(homogen,PolF5) # module de la solution homogène
H2 = [z2,mod2]


H,L,U = MySNF(X3)
b33 = L*b3
non_nul = [j for j in range(min(X3.nrows(),X3.ncols())) if H[j,j]!=0]
nul = [j for j in range(min(X3.nrows(),X3.ncols())) if H[j,j]==0]
sol_part = []
for k in non_nul:
    sol_part.append((b33[k]//H[k,k])*U[:,k])
z3 = sum(sol_part)
homogen = [vector(U[:,k]) for k in nul]
mod3 = span(homogen,PolF5) # module de la solution homogène
H3 = [z3,mod3]
# # Affichage des resultats

print("Une solution particuliere de X1*z1 = b1 est")
print(z1)
print("les solutions du systeme homogene sont engendres par")
print(H1)
print("Une solution particuliere de X2*z2 = b2 est")
print(z2)
print("les solutions du systeme homogene sont engendrees par")
print(H2)
print("Une solution particuliere du systeme 3 est")
print(z3)
print("les solutions du systeme homogene sont engendres par")
print(H3)



#reset()
print("""\
# ****************************************************************************
# STRUCTURE DU QUOTIENT
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
              [ -630,   735,   0,   735, -630],
              [ 1275, -1485, -15, -1470, 1275],
              [  630,  -630,   0,  -630,  630]])

# Code pour l'EXERCICE
H,L,U = MySNF(A)
L_inv = L.inverse()
N = []
coeff = []
for k in range(min(H.ncols(),H.nrows())):
    N.append(vector(H[k,k]*L_inv[:,k]))
    coeff.append(H[k,k])
print(N)
print(coeff)
N = span(N,ZZ)
M = ZZ^3/N
print(M)
reponse = " ici M = (Z/15Z)*l1 + (Z/105Z)*l2 + (Z/630Z)*l3 avec l1,l2,l3 les colonnes de L^-1 sont: \
l1 = (42,-85,-42); \
l2 = (1,-2,0); \
l3 = (0,0,1)"

# # Affichage des resultats

print("La structure de Z^3/N est")
print(reponse)

print("""\
# ****************************************************************************
# FACTEURS INVARIANTS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
              [ -630,   735,   0,   735, -630],
              [ 1275, -1485, -15, -1470, 1275],
              [  630,  -630,   0,  -630,  630]])


# Code pour l'EXERCICE
H,L,U = MySNF(A)
L_inv = L.inverse()
N = []
for k in range(min(H.ncols(),H.nrows())):
    N.append(vector(H[k,k]*L_inv[:,k]))
N = span(N,ZZ) # génération du Z-module
M = ZZ^3/N
print(M)
rang = M.ngens() - len(M.invariants())
fact_inv = M.invariants()
exponents = fact_inv
reponse = "Le module est de rang 0, car le nombre d'invariants est égal à la dimension de l'espace, et implique que M est de torsion. Il possède trois éléments cycliques d'ordre respectifs 15, 105 et 630, avec 15|105 et 105|630: ces ordres sont ses exposants"

# Affichage des resultats

print("Le rang de Z^3 / N est")
print(rang)
print("Les facteurs invariants sont")
print(fact_inv)
print("Exposants ?")
print(reponse)


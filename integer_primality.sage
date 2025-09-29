print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP8 : PRIMALITE DES ENTIERS                                                 #
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
# TEST DE RABIN-MILLER 
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n = 561

# Code pour l'EXERCICE

def testRM(n,a=None):
    if n==2:
        return True
    a = ZZ.random_element(2,n) if a is None else a
    v = (n-1).valuation(2)
    m = (n-1)//2^v
    if gcd(a,n)>1:
        return False
    b = mod(a^m,n)
    if b==1:
        return True
    for i in range(1,v+1):
        if b^2.mod(n)==1:
            g = gcd(b+1,n)
            if g==1 or g==n:
                return True
            else:
                return False
        b = b^2.mod(n)
    #return ZZ(n).is_prime() # A MODIFIER
    return False

# # Affichage des resultats

print("Test de la primalite de n=",n,"avec implementation de Rabin-Miller")
print(testRM(n))

print("""\
# ****************************************************************************
#  PERFORMANCES DE RABIN-MILLER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

nmin=10
nmax=500
nbtests = 100

# Code pour l'EXERCICE

rep2 = "Les nombres composés déclarés premiers avec une probabilité supérieure à 0.1 sont: "
rep3 = "La proportion de témoins de Rabin-Miller est de l'ordre de 3/4, soit il faudrait 1/4^k pour k le nombre de tests échoués. Nous devons alors avoir au moins 25 tests différents."

# # Affichage des resultats

bar_chart( [sum( [testRM(n)/nbtests for i in range(nbtests)]) if not ZZ(n).is_prime() else 0 for n in range(nmin,nmax)]).save('barchart.png')

faux_pos = []
for n in range(nmin,nmax):
    if not ZZ(n).is_prime():
        if sum([testRM(n)/nbtests for i in range(nbtests)])>0.1:
            faux_pos.append(n)

print(rep2)
print(faux_pos)
print("Ces nombres sont composés de deux facteurs de nombres premiers.")
print(rep3)
list_plot( [timeit( 'testRM(n)', number=20, repeat=3, seconds=true) for n in range(1001,1001+100000,100) ]).save('exec.png')


#reset()
print("""\
# ****************************************************************************
# TEST DE SOLOVAY-STRASSEN 
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n = 561

# Code pour l'EXERCICE

def testSS(n,a=None):
    if mod(n,2)==0:
        return False
    a = ZZ.random_element(2,n) if a is None else a
    if gcd(a,n)!=1:
        return False
    if jacobi_symbol(a,n)==mod(a^((n-1)//2),n):
        return True
    else:
        return False
    #return ZZ(n).is_prime() # A MODIFIER


nmin=10
nmax=500
nbtests = 100


rep3 = "Les faux positifs de probabilité supérieure à 0.1 sont:"
rep4 = "Cette fois-ci, le nombre de témoins est de l'ordre de la moitié des candidats, il faut alors au moins k=50 tests avec des a distincts."

# # Affichage des resultats

print("Test de la primalite de n=",n,"avec implementation de Solovay-Strassen")
print(testSS(n))

bar_chart( [sum( [testSS(n)/nbtests for i in range(nbtests)]) if not ZZ(n).is_prime() else 0 for n in range(nmin,nmax)]).save('barchartSS.png')

faux_pos = []
for n in range(nmin,nmax):
    if not ZZ(n).is_prime():
        if sum([testSS(n)/nbtests for i in range(nbtests)])>0.1:
            faux_pos.append(n)


print(rep3)
print(faux_pos)
print("Ces nombres sont semblables aux faux positifs de Rabin Miller")
print(rep4)


#reset()
print("""\
# ****************************************************************************
# COMPARAISON ENTRE LES TESTS DE R-M ET S-S 
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

nmax=150

# Code pour l'EXERCICE

Temoins = []
for n in range(2,nmax):
    if not ZZ(n).is_prime():
        for a in range(2,n):
            if not testRM(n,a) and testSS(n,a):
                Temoins.append((n,a))

# # Affichage des resultats

print("Liste d'entiers composés et de temoins exclusifs de Rabin-Miller")
print(Temoins)



#reset()
print("""\
# ****************************************************************************
# TEST DE LUCAS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

nmax=1000

#p = 3
#q = 1


u = [-1 for i in range(nmax)]
v = [-1 for i in range(nmax)]

# def init_uv(p,q):
#     u = [-1 for i in range(nmax)]
#     v = [-1 for i in range(nmax)]
#     u[0] = 0
#     v[0] = 2
#     u[1] = 1
#     v[1] = p


# def calcul_u(k,p,q):
#     if u[k]!=-1:
#         return u[k]
#     if mod(k,2)==0:
#         u[k] = (calcul_u(k//2,p,q) if u[k//2]==-1 else u[k//2])*(calcul_v(k//2,p,q) if v[k//2]==-1 else v[k//2])
#     else:
#         u[k] = (calcul_u(k//2+1,p,q) if u[k//2+1]==-1 else u[k//2+1])*(calcul_v(k//2,p,q) if v[k//2]==-1 else v[k//2])-q^(k//2)
#     return u[k]

# def calcul_v(k,p,q):
#     if v[k]!=-1:
#         return v[k]
#     if mod(k,2)==0:
#         v[k] = (calcul_v(k//2,p,q)^2 if v[k//2]==-1 else v[k//2]^2)-2*q^(k//2)
#     else:
#         v[k] = (calcul_v(k//2+1,p,q) if v[k//2+1]==-1 else v[k//2+1])*(calcul_v(k//2,p,q) if v[k//2]==-1 else v[k//2])-p*q^(k//2)
#     return v[k]


def calcul_u_v(k,p,q):
    if k == 0:
        return 0,2
    if k == 1:
        return 1,p
    if k % 2 == 0:
        j = k/2
        uj,vj = calcul_u_v(j,p,q)     # k = 2j
        return uj*vj, vj^2 - 2*q^j
    if k % 2 == 1:
        j = (k-1)/2
        uj,vj = calcul_u_v(j,p,q)     # k = 2j+1
        uj1,vj1 =calcul_u_v(j+1,p,q)
        return uj1*vj - q^j, vj1*vj - p*q^j

def choisir_p_q(n):
    for p in range(3,100):
        for q in range(1,100):
            delta = p^2 - 4*q
            if delta == 0:
                continue
            if gcd(n, 2*q*delta) != 1:
                continue
            if jacobi_symbol(delta, n) == -1:
                return p, q
    print("Pas de (p,q) trouvés convenables dans l'intervalle.")
    return None,None


def testL(n,p,q):
    delta = p^2-4*q
    g = gcd(n,2*q*delta)
    if 1<g and g<n:
        return False
    elif g==n:
        p,q = choisir_p_q(n)
        if p is None:
            return False
        #raise ValueError("mauvais choix de p et q")
    t = (n-jacobi_symbol(delta,n)).valuation(2)
    m = (n-jacobi_symbol(delta,n))//2^t
    um,_ = calcul_u_v(m,p,q)
    g = gcd(n,um)
    if 1<g and g<n:
        return False
    elif g==n:
        return True
    for s in range(t):
        _,vs = calcul_u_v(2^s*m,p,q)
        g = gcd(n,vs)
        if 1<g and g<n:
            return False
        elif g==n:
            return True
    return False
    #return ZZ(n).is_prime() # A MODIFIER

# # Affichage des resultats

for _ in range(10):
    n =  ZZ.random_element(2,1000)
    if not n.is_prime():
        print(n.is_prime()==testL(n,3,1),n)

#reset()
print("""\
# ****************************************************************************
# TEST DE BAILLIE, POMERANCE, SELFRIDGE ET WAGSTAFF
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice



# Code pour l'EXERCICE

def calcul_delta(n):
    k=0
    delta = 5
    while jacobi_symbol(delta,n)!=-1:
        k = k+1
        delta = (-1)^k*(2*k+5)
    return delta

def testBPSW(n):
    if n==2:
        return True
    if not testRM(n,2):
        return False
    delta = calcul_delta(n)
    p = 1
    v[1] = p
    q = (1-delta)//4
    if not testL(n,p,q):
        if ZZ(n).is_prime():
            print("L",n)
        return False
    return True
    #return ZZ(n).is_prime() # A MODIFIER

# # Affichage des resultats

test = [ZZ(n).is_prime()==testBPSW(n) for n in range(2,nmax+1)]
print(all(test))

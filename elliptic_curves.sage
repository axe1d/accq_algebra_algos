print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP13 : COURBES ELLIPTIQUES                                                  #
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
# ADDITION DANS UNE COURBE ELLIPTIQUE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 61
Fp = FiniteField(p)
E = EllipticCurve(Fp,[1,0])
P = E.random_point()
Q = E.random_point()

# Code pour l'EXERCICE

def addition(P,Q):
    E = P.curve()
    Fp = E.base_ring()
    a = E.a4()
    b = E.a6()
    if P.is_zero():
        return Q
    elif Q.is_zero():
        return P
    elif P[0]==Q[0] and P[1]==-Q[1]:
        return E(0)
    if P==Q:
        if P[1]==0:
            return E(0)
        lamb = (3*P[0]**2+a)/2*P[1]
        xr = lamb**2-2*P[0]
        yr = -P[1]+lamb*(P[0]-xr)
    else:
        lamb = (Q[1]-P[1])/(Q[0]-P[0])
        xr = lamb**2-P[0]-Q[0]
        yr = -P[1]+lamb*(P[0]-xr)
    return E(xr,yr) #FAIRE LE CALCUL DU POINT

# # Affichage des resultats

print("Addition implémentée de ",P,Q," en affine:",addition(P,Q))
print("Méthode SageMath:",P+Q)


#reset()
print("""\
# ****************************************************************************
# COURBE DE L'ANSSI
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE
# https://www.legifrance.gouv.fr/jorf/id/JORFTEXT000024668816
ANSSI = "Agence nationale de la sécurité des systèmes d’information"
p = Integer(int('F1FD178C0B3AD58F10126DE8CE42435B3961ADBCABC8CA6DE8FCF353D86E9C03',16))
a = Integer(int('F1FD178C0B3AD58F10126DE8CE42435B3961ADBCABC8CA6DE8FCF353D86E9C00',16))
b = Integer(int('EE353FCA5428A9300D4ABA754A44C00FDFEC0C9AE4B1A1803075ED967B7BB73F',16))
E = EllipticCurve(FiniteField(p),[a,b])

# # Affichage des resultats

print("ANSSI signifie :",ANSSI)
print("La courbe recommandée est")
print(E)
print("p est premier ?",p.is_prime())
print("L'ordre de E est:",E.order())
secu = int(floor(1/2*(E.order().log(2)+1))) # compter Oe...
print("Sécurité de E:",secu)

#reset()
print("""\
# ****************************************************************************
# COMPTAGE DE POINTS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 2003
Fp = FiniteField(p)
a = 1929
b = 1178
E = EllipticCurve(Fp,[a,b])

while true:
    d=Fp.random_element()
    if not d.is_square():
        break


# Code pour l'EXERCICE

def comptage(Fp,a,b):
    return p+1+sum([legendre_symbol(x0**3+a*x0+b,p) for x0 in Fp])


def count_curves(Fp):
    curves = []
    for a in Fp:
        for b in Fp:
            if -16*(4*a**3+27*b**2)!=0:
                curves.append((a,b))
    return curves


p=7
Fp7=FiniteField(p)

orders = [EllipticCurve(Fp7,[a,b]).order() for a,b in count_curves(Fp7)]
orders.sort()
unique_orders = list(set(orders))
y = []
for x in unique_orders:
    y.append(orders.count(x))
y = [0]*(unique_orders[0])+y
frequence = [4,1,4]
bar_chart(y).save('freq_ordre_F7.png')


while true:
    d=Fp7.random_element()
    if not d.is_square():
        break

def count_d(Fp,d):
    cards = []
    for a in Fp:
        for b in Fp:
            #cards.append(comptage(Fp,a,b)+comptage(Fp,d**2*a,d**3*b))
            if -16*(4*a**3+27*b**2)!=0:
                cards.append(EllipticCurve(Fp,[a,b]).count_points()+EllipticCurve(Fp,[d**2*a,d**3*b]).count_points())
    return cards

list_cards = count_d(Fp7,d)
print("On observe que toutes ces sommes sont égales à:",list_cards[0],"pour F7, au nombre de :",len(list_cards))


# # Affichage des resultats
p = 2003
Fp = FiniteField(p)
a = 1929
b = 1178
E = EllipticCurve(Fp,[a,b])
print("Implémentation bête comptage de E:",comptage(Fp,a,b))
print("Méthode Sage:",E.count_points())
print("Le diagramme est symétrique.")


#reset()
print("""\
# ****************************************************************************
# FACTORISATION ECM
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

class FoundFactor(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

n = Zmod(2020)
n=2020

# Code pour l'EXERCICE

def division(x,y):
    n = y.modulus()
    try:
        return Zmod(n)(x/y)
    except ZeroDivisionError:
        d = Integer(y.lift()).gcd(Integer(n))
        return FoundFactor(d)


def addition(P,Q):
    E = P.curve()
    Fp = E.base_ring()
    a = E.a4()
    b = E.a6()
    if P.is_zero():
        return Q
    elif Q.is_zero():
        return P
    elif P[0]==Q[0] and P[1]==-Q[1]:
        return E(0)
    if P==Q:
        if P[1]==0:
            return E(0)
        lamb = division((3*P[0]**2+a),2*P[1])
        if isinstance(lamb,FoundFactor):
            raise lamb
        xr = lamb**2-2*P[0]
        yr = -P[1]+lamb*(P[0]-xr)
    else:
        lamb = division((Q[1]-P[1]),(Q[0]-P[0]))
        if isinstance(lamb,FoundFactor):
            raise lamb
        xr = lamb**2-P[0]-Q[0]
        yr = -P[1]+lamb*(P[0]-xr)
    return E(xr,yr) #FAIRE LE CALCUL DU POINT


def multiplication(lamb,P):
    E = P.curve()
    if lamb==0:
        return E(0)
    elif lamb==1:
        return P
    elif lamb%2==0:
        return multiplication(lamb//2,addition(P,P))
    else:
        return addition(multiplication(lamb//2,addition(P,P)),P)


def ECM(n,B):
    # a = n.random_element()
    # x0 = n.random_element()
    # y0 = n.random_element()
    a  = Zmod(n)(randint(0, n-1))
    x0 = Zmod(n)(randint(0, n-1))
    y0 = Zmod(n)(randint(0, n-1))
    n = a.modulus()
    b = y0**2-x0**3-a*x0
    g = Integer(4*a**3+27*b**2).gcd(Integer(n))
    if 1<g and g<n:
        return Integer(g)
    elif g==n:
        raise False
    E = EllipticCurve(Zmod(n),[a,b])
    A = E(x0,y0)
    for p in prime_range(B+1):
        cnt = 1
        while p**(cnt+1)<=B:
            cnt+=1
        try:
            A = multiplication(p**cnt,A)
        except FoundFactor as f:
            return Integer(f.value)
    return False

# # Affichage des resultats

print(ECM(n,15))


#reset()
print("""\
# ****************************************************************************
# EXPONENTIATION RAPIDE ET ATTAQUE PAR CANAUX CACHES
# ****************************************************************************
""")

# NE PAS TRAITER



#reset()
print("""\
# ****************************************************************************
# COURBE D'EDWARDS
# ****************************************************************************
""")



# Code pour l'EXERCICE

reponse = "Ces courbes empêchent les attaques par canaux auxiliaires du fait de la similitude entre les opérations d'addition et de duplication."

# # Affichage des resultats

print(reponse)


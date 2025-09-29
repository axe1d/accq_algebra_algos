import math
from sage.plot.plot3d.parametric_plot3d import parametric_plot3d
from sage.plot.plot3d.implicit_plot3d import implicit_plot3d
from sage.plot.plot3d.list_plot3d import list_plot3d
#from sage.plot.plot3d.shapes import line3d
import numpy as np
np.complex = np.complex128  # Temporary fix for deprecated alias

print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP6 : BASES DE GROEBNER ET SYSTEMES POLYNOMIAUX MULTIVARIES                 #
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
#  FONCTIONS DE SAGEMATH
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f = 2*x^2*y+7*z^3

# Code pour l'EXERCICE

print(x<y^2)
print(f.lt())
print(f.lc())
print(f.lm())

reponse  ="Les méthodes lt, lc, lm correspondent respectivement au terme dominant, le coefficient dominant et le monôme dominant.\nLes autres ordres monomiaux sont: degrevlex, deglex, invlex, neglex, negdegrevlex, degdeglex, wdegrevlex, wdeglex, negwdegrevlex."

# # Affichage des resultats

print("\n$1/ ", reponse)

#reset()
print("""\
# ****************************************************************************
# DIVISION MULTIVARIEE
# ****************************************************************************
""")




# Donnees de l'enonce de l'exercice

MPol.<x,y> = PolynomialRing(QQ,2, order='lex')
f  = -x^7 + x^6*y + 2*x^5 - 2*x^4*y - 5*x^2 + 3*x*y^3 + 5*x*y + 11*y^3 + 10 
f1 = x*y^2+2*y^2
f2 = x^5+5

# Code pour l'EXERCICE

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    Q = [MPol(0)]*s
    r = 0
    p = f
    while p!=0:
        found = False
        for i in range(s):
            if p.lt()%F[i].lt()==0:
                Q[i] = Q[i]+(p.lt()//F[i].lt())
                p = p-(p.lt()//F[i].lt())*F[i]
                found = True
                break
        if not found:
            r = r+p.lt()
            p = p-p.lt()
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r

# # Affichage des resultats

print("\n$ ",  myDivision(f,[f1,f2]))

#reset()
print("""\
# ****************************************************************************
# BASE DE GROEBNER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f1 = x^2-y
f2 = x*y-z
f3 = z^4+x*y

# Code pour l'EXERCICE

def syzygie(g,h):
    return lcm(g.lt(),g.lt())//g.lt()*g-lcm(g.lt(),h.lt())//h.lt()*h

def myGroebner(F):
    G = Set(F)
    S = Set()
    cnt = 0
    while not S.is_empty() or cnt==0:
        S = Set()
        for i in range(len(F)):
            for j in range(len(F)):
                if i!=j:
                    r = myDivision(syzygie(G[i],G[j]),G)[1]
                    if r!=0:
                        S = S.union(Set([r]))
        cnt+=1
        G = G.union(S)
    return G

def myRedGroebner(F):
    G = myGroebner(F)
    G_l = list(G)
    cond = False
    while not cond:
        cond = True
        for i in range(len(list(G))):
            H = list(G.difference(Set([G_l[i]])))
            G_l[i] = myDivision(G_l[i].lt(),H)[1]
            if G_l[i]==0:
                cond = False
            G_l[i] = G_l[i]//G_l[i].lc()
    return Set(G_l)

# # Affichage des resultats

print("\n$1/ ",myGroebner([f1,f2,f3]))
print("\n$2/ ",myRedGroebner([f1,f2,f3]))




#reset()
print("""\
# ****************************************************************************
# APPARTENANCE A UN IDEAL
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f1 = x*y-y^2
f2 = x^3-z^2
I = Ideal([f1,f2])
f = -4*x^2*y^2*z^2 + y^6 + 3*z^5

# Code pour l'EXERCICE

test1 = f in I
test2 = myDivision(f,list(myRedGroebner([f1,f2])))[1]==0  # A ECRIRE VOUS-MEME


# # Affichage des resultats

print("\n$ Test de Sage ",test1)
print("\n$ Test de personnel ",test2)


#reset()
print("""\
# ****************************************************************************
# RESOLUTION D'UN SYSTEME
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


MPol.<x,y> = PolynomialRing(QQ,2, order='lex') # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
f = (y^2+6)*(x-1) - y*(x^2 + 1)
g = (x^2+6)*(y-1) - x*(y^2 + 1)


# Code pour l'EXERCICE

base = Ideal([f,g]).groebner_basis() # Vous pouvez utiliser la fonction adhoc de sage
          # pour calculer la base Groebner
# Le deuxième polynôme est le seul faisant intervenir des monômes en x et y

racines_y = base[2].univariate_polynomial().roots()
racines_y = [x[0] for x in racines_y]
racines = []
l1 = [(z,base[0].subs({y:z}).univariate_polynomial().roots()) for z in racines_y if base[0].subs({y:z}).univariate_polynomial() != 0]
l2 = [(z,base[1].subs({y:z}).univariate_polynomial().roots()) for z in racines_y if base[1].subs({y:z}).univariate_polynomial() != 0]

for (z,k) in l1:
    racines = racines+[(z,x[0]) for x in k]

for (z,k) in l2:
    racines = racines+[(z,x[0]) for x in k]

Gf = implicit_plot(f,(x,0,6),(y,0,6),color='red')
Gg = implicit_plot(g,(x,0,6),(y,0,6),color='blue')
Gp = point2d(racines,color='green')

# # Affichage des resultats

print("\n$1/  Une base de Groebner de [f,g] est", base)
print("\n$2/  Les valeurs de y sont", racines_y)
print("\n$4/  Les valeurs de (x,y) sont", racines)
print("\n$4/")
show(Gf+Gg+Gp)
(Gf+Gg+Gp).save('syst_resol.png')


#reset()
print("""\
# ****************************************************************************
# OPTIMISATION
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


MPol.<x,y,lamb> = PolynomialRing(QQ,3,order='lex') # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
f = x^2*y  - 2*x*y + y + 1
g = x^2 + y^2 - 1

# Code pour l'EXERCICE

syst = [2*x*y-2*y-lamb*(2*x),x^2-2*x+1-lamb*(2*y),g]
base = Ideal(syst).groebner_basis()
#racines_lamb = base[3].univariate_polynomial().roots()
#racines_lamb = [z[0] for z in racines_lamb if z[0]!=0]
racines_lamb = [math.sqrt(125/36),-math.sqrt(125/36)]
racines_x = []
racines_y = []
racines = []

lx = [base[1].subs({lamb:z}).univariate_polynomial().roots() for z in racines_lamb]
for k in lx:
    racines_x = racines_x+[z[0] for z in k]

ly = [base[2].subs({lamb:z}).univariate_polynomial().roots() for z in racines_lamb]
for k in ly:
    racines_y = racines_y+[z[0] for z in k]

for i in range(len(racines_x)):
    racines.append((racines_x[i],racines_y[i]))

# # Affichage des resultats


print("\n$1/  On doit resoudre le systeme", syst)
print("\n$2/  dont une base de Groebner est", base)
print("\n$4/  Les valeurs de (x,y) sont", racines)

# Paramétrisation tracé de g
theta = var('theta')
x = cos(theta)
y = sin(theta)
g_theta = x^2 * y - 2 * x * y + y + 1
plot_g = plot(g_theta, (theta, 0, 2*pi), color='blue', legend_label='g(θ)')
plot_g.show()
plot_g.save('opti.png')

#reset()
print("""\
# ****************************************************************************
# MANIPULATIONS ALGEBRIQUES
# ****************************************************************************
""")
print("""\
# ****************************************************************************
# Exercice 261
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice
MPol.<x,y,z,t> = PolynomialRing(QQ,4,order='neglex') # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
syst = [x-t^2,y-t^3,z-t^4]
base = Ideal(syst).groebner_basis()
# Equation y^2 = z*x

# Code pour l'EXERCICE
A.<x,y,z> = AffineSpace(QQ,3)
C = Curve([y^2-x*z])

x, y, z, t = var('x y z t')

# Implicit plot for y^2 = x*z
implicit_curve = implicit_plot3d(
    y^2 - x*z,
    (x,0,25), (y,-125,125), (z,0,625),
    color='red',
    opacity=0.4,
    plot_points=100
)

# Parametric plot for the curve
param_curve_rational = parametric_plot3d(
    (t^2, t^3, t^4),  # x = t^2, y = t^3, z = t^4
    (t, -5, 5),
    color='blue',
    thickness=2,
    legend_label='Parametric Curve'
)

combined_plot = implicit_curve + param_curve_rational

combined_plot.save('curve.png')

# # Affichage des resultats

print("\n$1/  On doit resoudre le systeme", syst)
print("\n$2/  dont une base de Groebner est", base)
print("\n$4/  Les valeurs de (x,y,z) sont données par l'équation y^2 = x*z\n")



print("""\
# ****************************************************************************
# Exercice 263
# ****************************************************************************
""")

MPol.<x,y,z,t> = PolynomialRing(QQ,4,order='lex') # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
syst = [x^2+y^2-1,z^2+t^2-1,z-2*x*y,t-y^2+x^2,t-2*y^2+1]
base = Ideal(syst).groebner_basis()

u,v = var('u v')
syst = [(f.subs({y:u-x,t:v-z})) for f in base]

MPol_n.<x,z,u,v> = PolynomialRing(QQ,4,order='lex') # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
syst = [MPol_n(f) for f in syst]
base = Ideal(syst).groebner_basis()

red = MPol_n(x^6).reduce(base)

print("\nOn obtient finalement sin^6(theta) =",red,"avec u = sin(theta)+cos(theta) et v = sin(2theta)+cos(2theta)")

#reset()
print("""\
# ****************************************************************************
# OVALES DE DESCARTES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice
MPol.<x,y,z,t> = PolynomialRing(QQ,4,order='lex') # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
# OM+2O'M = 3
# OM = z; O'M = t
# O = (0,0) ; O' = (1,0)
syst = [z+2*t-3,x^2+y^2-z^2,(x-1)^2+y^2-t^2]
base = Ideal(syst).groebner_basis()

# z = -2t+3
base = [f.subs({z:base[-1]}) for f in base[:-1]]
t = var('t')
solutions = solve(x-(3/2)*t^2+6*t-5==0,t)
print(solutions)
x,y = var('x y')
t = -1/3*(6*x+6)^(1/2)+2
f = y^2+2.25*t^4-18*t^3+47*t^2-48*t+16
Gf = implicit_plot(f,(x,-1,2),(y,-1,1),color='red')
show(Gf)
Gf.save('descartes.png')

# Code pour l'EXERCICE

eq = f.simplify_full()

# # Affichage des resultats

print("\n$1/  On doit resoudre le systeme", syst)
print("\n$2/  dont une base de Groebner est", base)
print("\n$ L'équation est ",eq)

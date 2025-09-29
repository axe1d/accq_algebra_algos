from sage.plot.plot3d.parametric_plot3d import parametric_plot3d

print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP7 : APPLICATIONS CHOISIES DES BASES DE GROEBNER                           #
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
# POINTS SINGULIERS 
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y> = PolynomialRing(QQ,2,order='lex')
f = 5*x^4 - 10*x^3 + 10*x^2*y^2 - 40*x^2*y + 40*x^2 - 10*x*y^2 + 40*x*y - 32*x + 5*y^4 - 40*y^3 + 115*y^2 - 136*y + 48

# Code pour l'EXERCICE
g1 = f.derivative(x)
g2 = f.derivative(y)
syst = [f,g1,g2]
b = Ideal(syst).groebner_basis()

y_sol = b[1].univariate_polynomial().roots()[0][0]
x_sol = b[0].subs({y:y_sol}).univariate_polynomial().roots()[0][0]
PtsSinguliers = [(x_sol,y_sol)] # A calculer

# # Affichage des resultats

print("La liste des points d'inflexions est :",PtsSinguliers)
G1 = implicit_plot(f,(x,-2,6),(y,-2,6),color='blue')
G2 = points(PtsSinguliers,color='red')
show(G1+G2)
#(G1+G2).save('points_singuliers.png')

#reset()
print("""\
# ****************************************************************************
#  VALUATIONS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

A.<x,y,z> = AffineSpace(QQ, 3)
QQ5.<rac5> = QuadraticField(5)
f1 = x^2+y^2-1
f2 = 5*x-(z-3)^2
Z = Curve([f1,f2],A)
n = z^2-6*z+5
d = x^2-x+y^2

# Code pour l'EXERCICE

Z1a = implicit_plot3d(
    x^2+y^2-1,
    (x,0,1), (y,-1,1), (z,3,5),
    color='red',
    opacity=0.4,
    plot_points=100
)

Z1b = implicit_plot3d(
    5*x-(z-3)^3,
    (x,0,1), (y,-1,1), (z,3,5),
    color='green',
    opacity=0.4,
    plot_points=100
)


def test_deg(A,R,t):
    for a in A:
        if all(t[i] >= a[i] for i in range(len(t))):
            return False
    return True


def dim_quo(I,n):
    I = Ideal(I).groebner_basis()
    A = Set([f.degrees() for f in I])
    R = Set([tuple([0 for _ in range(n)])])
    for i in range(n):
        S = R
        while not S.is_empty():
            t = S.an_element()
            t_list = list(t)
            S = S.difference(Set([t]))
            t_list[i] = t_list[i]+1
            while test_deg(A,R,t_list):
                if not tuple(t_list) in R:
                    R = R.union(Set([tuple(t_list)]))
                    t_list[i] = t_list[i]+1
    return R.cardinality()

zeros_n = n.univariate_polynomial().roots()
n1 = (z-zeros_n[0][0])^6
n2 = (z-zeros_n[1][0])^6
zeros = [((x,y,zeros_n[0][0]),dim_quo([n,f1,f2,n1],3)),((x,y,zeros_n[1][0]),dim_quo([n,f1,f2,n2],3))]

b = Ideal([d,f1,f2]).groebner_basis()
rac_x = b[2].univariate_polynomial().roots()[0][0]
poles_n = d.subs({x:rac_x}).univariate_polynomial().roots()
n1 = (y-poles_n[0][0])^6
poles = [((rac_x,poles_n[0][0],z),dim_quo([d,f1,f2,n1],3))]

var('theta')
x(theta) = theta
y(theta) = theta^2
z(theta) = theta^3

Z2 = parametric_plot3d(
    (sin(theta), cos(theta), 3+sqrt(5*sin(theta))),
    (theta, 0, pi),
    color='blue',
    thickness=5
)
#Z2.save('Z1.png')
#show(Z+Z2)
#(Z1a+Z1b).save('test.png')
#Z2.save('param_theta.png')
#(Z+Z2).save('cylindre2.png')

syst = Ideal(list(vector(f1.gradient()).cross_product(vector(f2.gradient()))))
base = syst.groebner_basis()
print("Base de Groebner du système :",base)
cond_lisse = not ((0,3,0) in Z)


reponse2="La courbe Z est lisse ? <=> le point (0,3,0) n'appartient pas à Z ? :"+str(cond_lisse)


# # Affichage des resultats

#parametric_plot3d((x,y,z),(theta,0,2)).save('result.png')


print(reponse2)
print("Confirmation de la lissité de la courbe : ", Z.is_smooth())

print("L'ensemble des zeros munis de leur multiplicité de h est",zeros)
print("L'ensemble des poles munis de leur multiplicité de h est",poles)



#reset()
print("""\
# ****************************************************************************
#  ENVELOPPE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice
MPol.<x,y,t> = PolynomialRing(QQ,3,order='invlex') # ORDRE DES VARIABLES A FIXER
f = (x-t)^2+(y+t^2-4*t)^2-4
g(x,y) = f.subs({t:1})
g2(x,y) = f.subs({t:2})
g3(x,y) = f.subs({t:3})

# Code pour l'EXERCICE
syst = [f,f.derivative(t)]
base = Ideal(syst).groebner_basis()
reponse1 = "La courbe (Ct) est un cercle dont le centre décrit une parabole quand t varie"

eq_enveloppe = base[-1]
f(x,y) = eq_enveloppe

l = y-x+5
syst = [eq_enveloppe,MPol(l)]
b = Ideal(syst).groebner_basis()

reponse4 = len(b[1].univariate_polynomial().roots())>0 # a-t-on un zéro commun ?

# # Affichage des resultats
print(reponse1)
print("L'enveloppe de la famille de courbe a pour equation",eq_enveloppe)
Z1 = implicit_plot(f,(x,-5,5),(y,-5,5),color='blue')
Z2 = implicit_plot(g,(x,-5,5),(y,-5,5),color='red')
Z3 = implicit_plot(g2,(x,-5,5),(y,-5,5),color='red')
Z4 = implicit_plot(g3,(x,-5,5),(y,-5,5),color='red')
show(Z1+Z2+Z3+Z4)
#(Z1+Z2+Z3+Z4).save('toit.png')


# Question 4


print("La bille expulsée tape-t-elle le toit ?", reponse4)

#reset()
print("""\
# ****************************************************************************
#  COLORATION DE GRAPHES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

G= Graph(12)
G.add_cycle(range(8))
G.add_edges([(i,i+4) for i in range(4) ])
G.add_edges([(8,5),(8,7),(9,0),(9,2),(10,1), (10,3),(11,4),(11,6) ])
G.add_edges([(8,9),(9,10),(10,11)])
G.show()
G.coloring()

MPol = PolynomialRing(QQ,12,'x',order = 'invlex')

K.<j> = CyclotomicField(3) # Pour coloration de graphe, travailler avec le groupe des racines 3-ième de l'unité
CPol = PolynomialRing(K,12,'x',order='invlex')

phi (v) = v^3-1
psi (u,v) = u^2+u*v+v^2
IG= Ideal(MPol, [phi(MPol.gen(v)) for v in G.vertices()] + [psi(MPol.gen(u),MPol.gen(v)) for (u,v) in G.edges(labels=false)])

# Code pour l'EXERCICE

reponse1 = 3
base = IG.groebner_basis()


# Coloration du graphe
color = [0 for _ in range(12)]

for i in range(11,-1,-1):
    d_subs = {MPol.gens()[j]:color[j] for j in range(0,11-i)}
    color[11-i] = base[i].subs(d_subs)
    if CPol(color[11-i]).degree()>0:
        color[11-i] = CPol(color[11-i]).univariate_polynomial().roots()[0][0]

reponse3 = "Le système est triangulaire avec 12 équations, résolvons le:"
reponse3 += "\nVoici la liste des valeurs de la coloration:" + str({MPol.gens()[j]:color[j] for j in range(12)})
reponse3 += "\nNous avons trois valeurs possibles, et la coloration est unique par invariance selon la valeur de départ donnée à x0"

# # Affichage des resultats

print("Il faut",reponse1,"couleurs pour colorer ce graphe")
print("Une base de Groebner de I(G,3) est", base)
print(reponse3)

#reset()
print("""\
# ****************************************************************************
#  PREUVE DE THEOREMES GEOMETRIQUES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y,u,v> = PolynomialRing(QQ,4)


# Code pour l'EXERCICE

eq_A = (u+1)*y - v*x
eq_B = (u-2)*y - v*(x-1)
eq_C = (2*u-1)*y - 2*v*(x-1/2)
eq_C = (2*u-1)*y - v*(2*x-1)


IA = Ideal(MPol,[eq_A])
IB = Ideal(MPol,[eq_B])
IC = Ideal(MPol,[eq_C])


IAB = IA+IB
IAetB = Ideal(MPol,[eq_A,eq_B])
reponse2  = "A-t-on IA + IB = I(A^B) ? " + str(IAB==IAetB)

reponse3 = "A-t-on IC inclus dans I(A^B) ? " + str(IAetB.intersection(IC)==IC)
reponse4 = "Nous venons de prouver la concurrence des médianes d'un triangle en un centroïde, par l'annulation du polynôme générateur de l'idéal de l'une sur l'intersection des deux autres idéaux."

# # Affichage des resultats

print("Les idéaux IA, IB et IC sont", IA, IB, IC)
print(reponse2)
print(reponse3)
print(reponse4)



#reset()
print("""\
# ****************************************************************************
#  PROGRAMMATION ENTIERE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<p,n,d,q> = PolynomialRing(QQ,4,order='degrevlex')
I=Ideal([p^5-n,p^10-d,p^25-q])

# Code pour l'EXERCICE

base = I.groebner_basis()
reponse2 = "L'intérêt de cette réduction modulo I est de décomposer dans l'anneau avec un ordre lexicographique, i.e p<n<d<q, afin de décomposer p^117 comme produits de facteurs en p^5 (nickel), p^10 (dime), p^25 (quarter) ainsi qu'en p^1 (penny). L'ordre implique la minimisation du nombre de pièces utilisées en privilégiant celles de plus forte valeur.\nOn obtient au final la décomposition finale après réduction: 4 quarters, 1 dime, 1 nickel et 2 pennys."

# # Affichage des resultats

print("Base de Groebner",base)
print((p^117).reduce(I))
print(reponse2)

#reset()
print("""\
# ****************************************************************************
#  SURFACE DE CLEBSCH
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x1,x2,x3> = PolynomialRing(QQ,3)
f=x1^3+x2^3+x3^3+1-(x1+x2+x3+1)^3

# Code pour l'EXERCICE
G1 = implicit_plot3d(f,(x1,-5,5),(x2,-5,5),(x3,-5,5))

MPold.<a,b,c,d,e,f> = PolynomialRing(QQ,6,order='lex')
T.<t> = PolynomialRing(MPold)
MPol.<x0,x1,x2,x3,x4> = PolynomialRing(T,5,order = 'lex')
F = x0^3+x1^3+x2^3+x3^3+x4^3
G = x0+x1+x2+x3+x4

x0_t = T(1)
x1_t = t
x2_t = a+d*t
x3_t = b+e*t
x4_t = c+f*t

ft = F.subs({x0:x0_t,x1:x1_t,x2:x2_t,x3:x3_t,x4:x4_t})
gt = G.subs({x0:x0_t,x1:x1_t,x2:x2_t,x3:x3_t,x4:x4_t})


ft = T(ft)
gt = T(gt)
coeff_ft = ft.coefficients()
coeff_gt = gt.coefficients()
I = Ideal(coeff_ft+coeff_gt)
base = I.groebner_basis()
variety = I.variety()


line_plots = []
for s in variety:
    a_val = s['a']
    b_val = s['b']
    c_val = s['c']
    d_val = s['d']
    e_val = s['e']
    f_val = s['f']

    L = parametric_plot3d([t, a_val + d_val*t, b_val + e_val*t],
                           (t, -5, 5), color='red', thickness=5)
    line_plots.append(L)

# # Affichage des resultats
#G1.save('clebsch_3.png')
#(G1+sum(line_plots)).save('droites_clebsch.png')
print("La surface de R³ est un espace affine pouvant être relevé en une partie de la surface C de Clebsch dans R⁴.\nEn effet, en introduisant X0, x1=X1/X0, x2=X2/X0, x3 = X3/X0, on obtient par projection l'équation X1³+X2³+X3³+X0³ - (X1+X2+X3+X0)³.")


print("Base de Groebner obtenue en {a,b,c,d,e,f}:",base)
print("Variété obtenue:",variety)
print("Nous obtenons au total",len(variety),"droites réelles contenues dans la surface de Clebsch")

print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP3 : RESEAUX EUCLIDIENS                                                    #
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
import os
from tp2_daboust import *

#reset()
print("""\
# ****************************************************************************
# BASE D'UN RESEAU
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
        [ 2, -3,  4, 0],
        [ 4,  4, 53, -5]])

# Code pour l'EXERCICE
# On utilise la forme normale d'Hermite car solution noyau de matrice

H,U = SageHNF(A)
L = matrix(ZZ,[vector(U[:-1,k]) for k in range(H.ncols()) if vector(H[:,k]) == vector([0 for _ in range(H.nrows())])]).transpose() # on ne prend pas la dernière coordonnée car correspond aux multiples de 5

# # Affichage des resultats

print("\n$ Le réseau a pour base")
print(L)



print("""\
# ****************************************************************************
# APPLICATIONS NUMERIQUES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice
R = RealField(200)
n1 = round(arctan(1),50)
n2 = round(arctan(1/5),50)
n3 = round(arctan(1/239),50)

r=-2.5468182768840820791359975088097915

Pol.<x>=PolynomialRing(ZZ)

M = 10^8

# Code pour l'EXERCICE
A = matrix(ZZ,[
   [1,0,0],
   [0,1,0],
   [0,0,1],
   [round(M*n1),round(M*n2),round(M*n3)]])
alpha = transpose(A).LLL()[0]

A = matrix(ZZ,[
   [1,0,0,0],
   [0,1,0,0],
   [0,0,1,0],
   [0,0,0,1],
   [round(M),round(M*r),round(M*r^2),round(M*r^3)]])
coeff = transpose(A).LLL()[0]

p = Pol(coeff[0] + coeff[1]*x + coeff[2]*x^2 + coeff[3]*x^3)

# # Affichage des resultats

print("\n$ La relation de Machin est alpha1*n1+alpha2*n2+alpha3*n3=0 avec")
for i in range(3):
   print("alpha",i+1,"=",alpha[i])

print("\n$ Un polynome minimal plausible est")
print(p)
print("dont les racines sont")
print(p.roots(ring=RR,multiplicities=false))




print("""\
# ****************************************************************************
# ALGORITHME LLL
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

B = matrix(ZZ,[[  9, -25],
               [ -5,  14],
               [  4, -11]])

# Code pour l'EXERCICE

def myLLL(M):
   B = copy(M)
   n = B.ncols()
   repeter = True
   while repeter:
      Betoile, mu = B.transpose().gram_schmidt(orthonormal=false)
      Betoile = Betoile.transpose()
      mu = mu.transpose()
      # ECRIRE L'ETAPE DE REDUCTION
      for i in range(1,n):
         for j in range(i-1,-1,-1):
            B.add_multiple_of_column(i,j,-round(mu[j,i]))
            mu.add_multiple_of_column(i,j,-round(mu[j,i]))
      assert(B==Betoile*mu)
      assert(all(mu[i,j]<=1/2 for i in range(n) for j in range(i+1,n)))
      assert(all(mu[i,j]>=-1/2 for i in range(n) for j in range(i+1,n)))

      # ECRIRE LE TEST D'ECHANGE
      echange = False
      i0 = 0
      while i0<n-1 and not echange:
         if vector(B[:,i0]).norm()^2 > 2*vector(B[:,i0+1]).norm()^2:
            echange = True
            B.swap_columns(i0,i0+1)
         i0 += 1

      if not echange:
         repeter = False
   return B


def myLLL_improved(M):
   B = copy(M)
   n = B.ncols()
   repeter = True
   Betoile, mu = B.transpose().gram_schmidt(orthonormal=false)
   Betoile = Betoile.transpose()
   mu = mu.transpose()
   while repeter:
      # ECRIRE L'ETAPE DE REDUCTION
      for i in range(1,n):
         for j in range(i-1,-1,-1):
            B.add_multiple_of_column(i,j,-round(mu[j,i]))
            mu.add_multiple_of_column(i,j,-round(mu[j,i]))
 
      assert(B==Betoile*mu)
      assert(all(mu[i,j]<=1/2 for i in range(n) for j in range(i+1,n)))
      assert(all(mu[i,j]>=-1/2 for i in range(n) for j in range(i+1,n)))

      # ECRIRE LE TEST D'ECHANGE
      echange = False
      i0 = 0
      while i0<n-1 and not echange:
         if vector(B[:,i0]).norm()^2 > 2*vector(B[:,i0+1]).norm()^2:
            echange = True
            B.swap_columns(i0,i0+1)
            s = Betoile[:,i0]
            t = Betoile[:,i0+1] + mu[i0,i0+1]*Betoile[:,i0]
            Betoile[:,i0] = t
            mu[i0,i0+1] = vector(s).dot_product(vector(Betoile[:,i0]))/(vector(Betoile[:,i0]).norm()^2)
            Betoile[:,i0+1] = s - mu[i0,i0+1]*Betoile[:,i0]
            for k in range(i0+2,n):
               mu[i0,k] = vector(B[:,k]).dot_product(vector(Betoile[:,i0]))/(vector(Betoile[:,i0]).norm()^2)
               mu[i0+1,k] = vector(B[:,k]).dot_product(vector(Betoile[:,i0+1]))/(vector(Betoile[:,i0+1]).norm()^2)
         i0 += 1

      if not echange:
         repeter = False
   return B


# # Affichage des resultats

print("\n$ Une base LLL de B est")
print(myLLL(B))
print("\n$ Une base LLL de B sans recalculer Gram-Schmidt est")
print(myLLL_improved(B))


print("""\
# ****************************************************************************
# RESEAUX CLASSIQUE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

n = 8
e = vector([1/2]*8)

# Code pour l'EXERCICE
def gram_matrix(A):
   G = matrix(ZZ,A.ncols(),A.ncols())
   for i in range(A.ncols()):
      for j in range(A.ncols()):
         G[i,j] = vector(A[:,i]).dot_product(vector(A[:,j]))
   return G



A = matrix(ZZ,[1 for i in range(n+1)])
H,U = SageHNF(A)
An = matrix(ZZ,[vector(U[:,k]) for k in range(H.ncols()) if vector(H[:,k]) == vector([0 for _ in range(H.nrows())])]).transpose()
an = gram_matrix(An).det() #determinant de An

D = matrix(ZZ,[1 for i in range(n)]+[-2])
H,U = SageHNF(D)
Dn = matrix(ZZ,[vector(U[:-1,k]) for k in range(H.ncols()) if vector(H[:,k]) == vector([0 for _ in range(H.nrows())])]).transpose() # on ne prend pas la dernière coordonnée car sert de modulo 2
dn = gram_matrix(Dn).det()

E = matrix(QQ,[1/2 for _ in range(8)])
E8 = Dn.transpose().stack(E).transpose() # vecteurs dépendants: on enlève la deuxième colonne
E8 = matrix(QQ,[vector(E8[:,k]) for k in range(E8.ncols()) if k!=1]).transpose()
e8 = sqrt(gram_matrix(E8).det())

# Question 5
# conditions respectées pour vecteurs de la base => vrai pour tout le réseau car stabilité de la parité
test = []
for i in range(An.ncols()):
   for j in range(An.ncols()):
      if i==j:
         test.append(vector(An[:,i]).dot_product(vector(An[:,i]))%2==0)
      else:
         test.append(vector(An[:,i]).dot_product(vector(An[:,j]))%1==0)

for i in range(Dn.ncols()):
   for j in range(Dn.ncols()):
      if i==j:
         test.append(vector(Dn[:,i]).dot_product(vector(Dn[:,i]))%2==0)
      else:
         test.append(vector(Dn[:,i]).dot_product(vector(Dn[:,j]))%1==0)

for i in range(E8.ncols()):
   for j in range(E8.ncols()):
      if i==j:
         test.append(vector(E8[:,i]).dot_product(vector(E8[:,i]))%2==0)
      else:
         test.append(vector(E8[:,i]).dot_product(vector(E8[:,j]))%1==0)

test = all(test)
print("\n$ Question 5: Réseaux An, Dn et E8 pairs: ", test)

reponse6 = "Si présence du dernier vecteur dans la décomposition alors toutes demi-entières car coefficients sont Z-module, sinon toutes entières"
reponse7_an = len([vector(An[:,i]) for i in range(An.ncols()) if vector(An[:,i]).dot_product(vector(An[:,i])) == 2])
reponse7_an *= reponse7_an-1 + 2 # on peut pour chaque vecteur soustraire un autre de la base pour annuler la première composante et en avoir une nouvelle à -1 pour obtenir 2 en norme, ou alors juste prendre le vecteur ou son opposé
reponse7_dn = len([vector(Dn[:,i]) for i in range(Dn.ncols()) if vector(Dn[:,i]).dot_product(vector(Dn[:,i])) == 2])
reponse7_dn = reponse7_dn*(reponse7_dn-1+2)*2 # comme avant mais on peut aussi ajouter le premier vecteur à n'importe lequel des autres vecteurs et ajouter au lieu de soustraire le vecteur obtenu à celui considéré: nouvelle composante à 1 maintenant
reponse7_e8 = len([vector(E8[:,i]) for i in range(E8.ncols()) if vector(E8[:,i]).dot_product(vector(E8[:,i])) == 2])
reponse7_e8 = 2 + 2*reponse7_e8 + 2*(reponse7_e8+1)*(reponse7_e8) + 2*reponse7_e8*(reponse7_e8+1)
# plus ou moins le dernier vecteur, puis chaque vecteur minimal de la base, puis combinaisons avec coordonnées non entières (on compense la première coordonnée avec le premier vecteur pour ramener à +-1/2), puis pour coordonnées entières

reponse8 = "autant que de vecteurs de norme 2"

# # Affichage des resultats

print("\n$ Une base de An est")
print(An, "de déterminant",an)

print("\n$ Une base de Dn est")
print(Dn, "de déterminant",dn)

print("\n$ Une base de E8 est")
print(E8, "de déterminant",e8)

print("\n$ Vecteurs respectivement de An, Dn et E8 de norme carrée égale à 2: ",reponse7_an,reponse7_dn,reponse7_e8)


print("""\
# ****************************************************************************
# DENSITES OPTIMALES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE
def compute_density(n,A):
   H,U = SageHNF(A)
   A2 = matrix(ZZ,[vector(U[:,k]) for k in range(H.ncols()) if vector(H[:,k]) == vector([0 for _ in range(H.nrows())])]).transpose()
   det = gram_matrix(A2).det()
   shortest_vec = (myLLL_improved(A2)[:,0]).norm() # plus petit vecteur en norme de L
   density = (pi()^(n/2)/(gamma(n/2 + 1) * 2^n))*(shortest_vec^(n/2))/sqrt(det)
   return density

n = 2
A = matrix(ZZ,[1 for i in range(n+1)])
a2 = compute_density(n,A)

n = 3
A = matrix(ZZ,[1 for i in range(n+1)])
a3 = compute_density(n,A)

n = 4
D = matrix(ZZ,[1 for i in range(n)]+[-2])
d4 = compute_density(n,D)

n = 5
D = matrix(ZZ,[1 for i in range(n)]+[-2])
d5 = compute_density(n,D)

n = 8
e8 = gram_matrix(E8).det()
e8 = pi()^(n/2)/(gamma(n/2+1)*2^n)*(vector(myLLL_improved(E8)[:,0]).norm())^(n/2) / sqrt(e8)

# # Affichage des resultats

print("\n$ La densité de A2 est",a2)
print("\n$ La densité de A3 est",a3)
print("\n$ La densité de D4 est",d4)
print("\n$ La densité de D5 est",d5)
print("\n$ La densité de E8 est",e8)





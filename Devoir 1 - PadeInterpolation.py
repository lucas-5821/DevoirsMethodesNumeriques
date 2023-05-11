#
# PYTHON for DUMMIES 22-23
# Problème 1
#
# Interpolation de Padé
# Interpolation avec des quotients de polynômes
#
# Vincent Legat 
# Nathan Coppin
#
# -------------------------------------------------------------------------
# 

from numpy import *
from numpy.linalg import solve

# ============================================================
# FONCTIONS A MODIFIER [begin]
#
#

def padeInterpolationCompute(X,U):
    A = []
    u = 0
    for i in range (len(X)) :
        m = [1]
        for j in range(1,((len(X)-1)//2)+1) :
            m.append(X[u]**j)
        for k in range(1,((len(X)-1)//2)+1) :
            m.append(-U[u]*X[u]**k)
        u += 1
        A.append(array(m))
    a = solve(array(A),U)
    return a
 
def padeEval(a,x) :
    n = 0
    d = 1
    for i in range((len(a)//2)+1) :
        n += a[i]*x**i
    k = 1
    for j in range((len(a)//2)+1,len(a)) :
        d += a[j]*x**k
        k+= 1 
    uh = n/d
    return uh
  
  
    
#
# FONCTIONS A MODIFIER [end]
# ============================================================

#
# -1- Test de la fonction interpolation
#     On considère un jeu des 3 fonctions u(x)
#


def main() :

  n = 2
  u =  lambda x : cos(x)   
  X = linspace(-2,2,(2*n+1)) 
  U = u(X)

#
# -1- Calcul des coefficients de l'interpolation
#
  print("==== Computing the Padé approximation :-)") 
  a = padeInterpolationCompute(X,U)
  print(" a = ",list(a))

  
#
# -2- Evaluation l'interpolation de Padé
#     et de l'interpolation polynomiale de Lagrange
#
  x = linspace(-5,5,100)
  upade = padeEval(a,x)
  uh = polyval(polyfit(X,U,len(X)-1),x)

#
# -3- Et un joli plot :-)
#

  from matplotlib import pyplot as plt
  plt.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'silver'

  plt.figure('Interpolation de Padé n = %d ' % n)
  plt.plot(x,u(x),'-b',label='cosinus')
  plt.plot(x,uh,'-g',label='Lagrangian interpolation') 
  plt.plot(x,upade,'-r',label='Padé interpolation')
  plt.plot(X,U,'ob')
  plt.xlim((-5,5)); plt.ylim((-1.2,1.2))
  plt.legend(loc='upper right') 
  plt.show()


main()


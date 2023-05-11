#
# PYTHON for DUMMIES 22-23
# Problème 3
#
# Splines cubiques périodiques
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 
 
 
from numpy import *
from numpy.linalg import solve
 
 
# ============================================================
# FONCTIONS A MODIFIER [begin]
#
#
 
def spline(x, h, U):
 n = size(U)
 X = arange(0, n+1)*h
 i = searchsorted(X[1:],x) 
 
 A = identity((n))
 for k in range(n):
     A[k][k] = 4  
     if k != n-1 :
         A[k][k+1] = 1 ; A[k][k-1] = 1
 A[k][k-1] = 1 ; A[k][0] = 1 ; A[0][0] = 4

 b = U[range(-1, n-1)] - 2*U[0:n] + U[append(arange(1, n), [0])]
      
 u = solve(((h**2)/6)*A, b) 
 
 U = append(U, U[0])
 u = append(u, u[0])

 n = n+1
 p = u[0:n-1]/(6*h)
 q = (U[0:n-1]/h - u[0:n-1]*h/6)
 pmoins1 = u[1:n]/(6*h)
 qmoins1 = (U[1:n]/h - u[1:n]*h/6)

 return pmoins1[i]*(x - X[i])**3 + p[i]*(X[i+1] - x)**3 +  qmoins1[i]*(x - X[i]) + q[i]*(X[i+1] - x) 
 
#
#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
# -1- Interpolation d'un cercle :-)     
#
 
def main() :
 
  from matplotlib import pyplot as plt
  plt.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'lavender'
 
  n = 4;
  h = 3*pi/(2*(n+1));
  T = arange(0,3*pi/2,h)
  X = cos(T); Y = sin(T)
 
  fig = plt.figure("Splines cubiques et cercle :-)")
  plt.plot(X,Y,'.r',markersize=10)
  t = linspace(0,2*pi,100)
  plt.plot(cos(t),sin(t),'--r')
 
  t = linspace(0,3*pi/2,100)
  plt.plot(spline(t,h,X),spline(t,h,Y),'-b')
  plt.axis("equal"); plt.axis("off")
  plt.show()
 
if __name__ == '__main__': 
  main()
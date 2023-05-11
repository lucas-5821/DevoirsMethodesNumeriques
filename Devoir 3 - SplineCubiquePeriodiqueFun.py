#
# PYTHON for DUMMIES 22-23
# Problème 3
#
# Script de test un peu plus rigolo
# Pour introduite un point : faire un clic sur la figure
# Un double clic permet d'obtenir le calcul de la courbe d'Hermite
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 

import matplotlib
from matplotlib import pyplot as plt
from numpy import *
from splineTest import spline
from numpy.linalg import solve

# ====================== callback pour les événements avec la souris ======
#
#  Observer la gestion distincte du clic simple et double :-)
#  Apres un evenement, on redessine la figure avec draw()
#
def spline(x, h, U):
 n = size(U)
 X = arange(0, n+1)*h
 i = zeros(len(x),dtype=int)
 for j in range(1,n):
     i[X[j]<=x] = j
 
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
 

def mouse(event):
  global X,Y,n
  if (event.dblclick):
    t  = arange(0,n+0.001,0.001)
    x  = spline(t,1.0,X)
    y  = spline(t,1.0,Y)
    plt.plot(x,y,'-b')
    X,Y = [],[]; n = 0
  else :    
    x = event.xdata 
    y = event.ydata
    if (x != None and y != None) :
      n = n + 1
      X = append(X,[x])
      Y = append(Y,[y])
      print("New data : " + str(x) + "," + str(y))
      plt.plot([x],[y],'.r',markersize=10)
  fig.canvas.draw()


# ============================= mainProgram ===============================
 

matplotlib.rcParams['toolbar'] = 'None'
matplotlib.rcParams['lines.linewidth'] = 1
plt.rcParams['figure.facecolor'] = 'lavender'

X,Y = [],[]; n = 0   
fig = plt.figure("Cubic spline interpolation")
fig.canvas.mpl_connect('button_press_event',mouse)
plt.ylim((0,1)); plt.xlim((0,1.3)); plt.axis("off")

plt.show()
 






# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 22-23
# Problème 6
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

from numpy import *


def lorenz(Xstart,Xend,Ustart,n): 
  
  f = lambda x, u: array([10 * u[1] - 10 * u[0] ,28 * u[0] - u[1] - u[0] * u[2] ,u[0] * u[1] - 8 * u[2] / 3  ])
    
  X = linspace(Xstart,Xend,n+1) 
  U = zeros((n+1,3)) 
  U[0,:] = Ustart 
  h = (Xend - Xstart)/n 
  
  for i in range(n):  
      K1 = f(X[i]    ,U[i,:]       )
      K2 = f(X[i]+h/2,U[i,:]+K1*h/2)
      K3 = f(X[i]+h/2,U[i,:]+K2*h/2)
      K4 = f(X[i]+h  ,U[i,:]+K3*h  )
      U[i+1,:] = U[i,:] + h*(K1+2*K2+2*K3+K4)/6     
    
  return X,U


def main():
  
# ------------------------------------------------------------------------------------ 
#
# Script de test
#
#
# ------------------------------------------------------------------------------------



  from matplotlib import pyplot as plt
  plt.rcParams['toolbar'] = 'None'
  plt.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'lavender'
  plt.rcParams['axes.facecolor'] = 'lavender'


  plt.figure("Lorenz Equations")
  Xstart = 0; Xend = 100;
  Ustart = [0,1,0]
  n = 10000;

  X,U = lorenz(Xstart,Xend,Ustart,n)
  plt.plot(U[:,0],U[:,2],'-r',linewidth=0.5)
  plt.show()



main()  

# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 22-23
# Problème 10
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

from numpy import *
from timeit import default_timer as timer


from scipy.sparse.linalg import spsolve 
import scipy.sparse as sparse

def poissonSolve(nCut):
  n = 2*nCut + 1; m = n*n; h = 2/(n-1)  
  A = sparse.dok_matrix(sparse.eye(m),dtype=float64)
  B = zeros(m)  
  for i in range(1,n-1):
      for j in range(1,n-1):
         if (i > nCut or j < nCut): 
            index = i + j*n
            A[index,index] = 4
            A[index,index-1] = -1
            A[index,index+1] = -1
            A[index,index+n] = -1
            A[index,index-n] = -1
            B[index] = 1            
  A = A / (h*h)  
  return spsolve(A.tocsr(),B).reshape(n,n)


  

def main():
  n = 20
  tic = timer()
  U = poissonSolve(n) 
  print("      Elapsed time : %f seconds" % (timer() - tic))
  print(" ==== Maximum value of U : %10.8f " % amax(U))
  print(" ==== Minimum value of U : %10.8f " % amin(U))


  import matplotlib.pyplot as plt
  import matplotlib
  matplotlib.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'lavender'
  plt.rcParams['axes.facecolor'] = 'lavender'
  myColorMap = matplotlib.cm.jet


  X = linspace(-1,1,2*n+1); U = abs(U)
  plt.figure("Python as Matlab clone...")
  plt.contourf(X,X,U,10,cmap=myColorMap)
  plt.contour(X,X,U,10,colors='k',linewidths=1)
  plt.hlines(X,X.min(),X.max(),color='white',linewidths=0.5)
  plt.vlines(X,X.min(),X.max(),color='white',linewidths=0.5)
  plt.axis("off"); plt.axis("equal")
  plt.show()
        

main()      
    


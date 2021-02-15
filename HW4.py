import numpy as np
import scipy as sc
#Question 1





#Question 2
def arnoldi(A,k,b):
    Q = np.zeros((len(b),k))
    H = np.zeros((k+1,k))
    Q[:,0] = b/np.linalg.norm(b)
    for j in range(k):
        Q[:,j+1] = np.matmul(A,Q[:,j])
        for i in range(j+1):
            H[i,j] = np.dot(Q[:,i],Q[:,j+1])
            Q[:,j+1] = Q[:,j+1] - H[i,j]*Q[:,i]
        H[j+1,j] = np.linalg.norm(Q[:,j+1])
        Q[:,j+1] /= H[j+1,j]
    return H[:-1,:], Q
    

def mygmres(l,b,x0,n,M,A):
    Q = np.zeros((len(b),l+1))
    H = np.zeros((l+1,l))
    r0 = b - np.matmul(A,x0)
    for i in range(l):
        #Should this be here?
        #r0 = b - np.matmul(A,x0)
        e1 = np.zeros(l)
        e1[0] = 1
        H,Q = arnoldi(A,l,b)
        alpha = e1 * np.linalg.norm(r0)
        print("H = ", H)
        print("\nalpha = ", alpha, "\n")
        y, res, rnk, s = sc.linalg.lstsq(H,alpha)
    sol = np.matmul(Q[:,:i+1],y) + x0
    return sol

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import gmres
A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)

b = np.array([2, 4, -1])
x, exitCode = gmres(A, b)
print(x, 'Scipy rsolution')
A = A.toarray()
x = mygmres(23,b,x+1.5,3,[0],A)
print(x)
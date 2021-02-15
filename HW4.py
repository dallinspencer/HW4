import numpy as np

#Question 1





#Question 2
def arnoldi(A,k,b):
    Q = np.zeros((len(b),k))
    H = np.zeros((k+1,k))
    Q[:,0] = b/np.linalg.norm(b)
    for j in range(k-1):
        Q[:,j+1] = np.matmul(A,Q[:,j])
        for i in range(j+1):
            H[i,j] = np.dot(Q[:,i],Q[:,j+1])
            Q[:,j+1] = Q[:,j+1] - H[i,j]*Q[:,i]
        H[j+1,j] = np.linalg.norm(Q[:,j+1])
        Q[:,j+1] /= H[j+1,j]
    return H[:-1,:], Q
    

def mygmres(l,b,x0,n,M,A):
    Q = np.zeros((len(b),n+1))
    H = np.zeros((n+1,n))
    Q[:,0] = b/np.linalg.norm(b)
    x = b
    for i in range(l):
        for j in range(n):
            H, Q = arnoldi(A,n,x)
            e1 = np.zeros(n)
            e1[0]=1
            y, res, rnk, s = np.linalg.lstsq(H,np.linalg.norm(b)*e1)
            print("Q[:,j+1]=",Q[:,j])
            print("y=",y)
            x = np.matmul(Q[:,j+1],y)
    return x                     

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import gmres
A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)

b = np.array([2, 4, -1])
x, exitCode = gmres(A, b)
print(x, 'Scipy rsolution')
A = A.toarray()
x = mygmres(500,b,x+1.5,3,[0],A)
print(x)
import numpy as np

#Question 1





#Question 2
def arnoldi(A,n,x0):
    q = np.zeros((n,n))
    q[0,:] = x0/np.linalg.norm(x0)
    h = np.zeros((n,n))
    for i in range(0,n-1):
        v = np.matmul(A,q[i][:])
        
        for j in range(0,i-1):
            h[j][i] = np.dot(q,v)
            v = v-h[j][i]*q
        h[i+1][i] = np.linalg.norm(v)
        q[i+1][:] = v/h[i+1][i]
    return h,q
    

def mygmres(l,b,x0,n,M,A):
    Q = np.zeros((n,n))
    H = np.zeros((n,n))
    for i in range(l):
        x = x0
        for j in range(n):
            H, Q = arnoldi(A,n,x)
            H = np.transpose(H)
            e1 = np.zeros(n)
            e1[0]=1
            y, res, rnk, s = np.linalg.lstsq(H,np.linalg.norm(b)*e1)
            print("Q:",Q)
            print("y:",y)
            print("H:",H)
            print("Norm(b):",np.linalg.norm(b))
            x = np.matmul(Q,y)
    return x                     

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import gmres
A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)

b = np.array([2, 4, -1], dtype=float)
x, exitCode = gmres(A, b)
print(x, 'Scipy rsolution')
A = A.toarray()
x = mygmres(10,b,x+1.5,3,[0],A)
print(x)
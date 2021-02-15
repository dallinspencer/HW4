
A=randi([10 100], 100, 100);
n = 100;
x0 = A(3,:)';
M = eye(100);

if det(A) ~=0
    b  = randi([10 100],100, 1);
    maxit = 100;
    tol=1e-1;
    [sol,xs,ys,Vs,Hs] = gmres_matlab(A,b,maxit,x0, M, n);
end
diff = b - A*sol


A=randi([10 100], 100, 100);
n = 100;
A=randi([10 100], n, n);
x0 = A(3,:)';
M = eye(100);

if det(A) ~=0
    b  = randi([10 100],n, 1);
    maxit = 90;
    tol=1e-1;
    [sol,xs,ys,Vs,Hs] = gmres_matlab(A,b,maxit,x0, M, n);
end
diff = b - A*sol

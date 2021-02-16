% b = [2;4;-1];
% A = [[3,2,0];[1,-1,0];[0,5,1]];
% l = 25;
% n = 3;
% x0 = [2.5;-1.5;10.5];
% M = eye(3);
% 
% [xs,ys,Vs,Hs]=gmres_matlab(A,b,l,x0,M,n);
% xs
% ys
n = 100;
A=randi([10 100], n, n);
x0 = A(3,:)';
M = eye(100);
if det(A) ~=0
    b  = randi([10 100],n, 1);
    maxit = 90;
    tol=1e-1;
    %[xs,ys,Vs,Hs] = gmres_matlab(A,b,maxit,x0, M, n) ;
    [x, error, iter, flag] = gmres_try2(A,x0,b,M,maxit)
end

testl2
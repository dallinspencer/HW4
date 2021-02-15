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
A=randi([10 100], 100, 100);
n = 100;
x0 = A(3,:)';
M = eye(100);
if det(A) ~=0
    b  = randi([10 100],100, 1);
    maxit = 100;
    tol=1e-1;
    [xs,ys,~,~] = gmres_matlab(A,b,maxit,x0, M, n) ;
end
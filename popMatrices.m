function [A,b] = popMatrices(n, gamma)
    h = 1/n;    
    A = diag((2/h)*ones(1,n)) + diag((-1/h-gamma/2)*ones(1,n-1),-1) + diag((-1/h+gamma/2)*ones(1,n-1),1);
    b = zeros(n,1) + h;
    b(1) = h/2;
    b(end) = h/2;       
%gamma = 1
figure(1)
count = 1;
errMat = [];%[n l error]
sgtitle('FEM Solutions')
for n = [16,32,64,128]
    l = 2;
    error = 1;
    [A,b] = popMatrices(n,1);  
    x0 = ones(n,1);
    M = eye(n);
    while error > 10^(-6)
        [sol,xs,ys,Vs,Hs] = gmres_matlab(A,b,l,x0, M, n);
        res = b - A*sol;
        error = norm(res)/n;
        l = l * 2;
    end
    errMat(end+1,:) = [n l error];
    
    x = linspace(0,1,n);
    subplot(2,2,count)
    plot(x,sol)
    title(sprintf('Solution with n = %d and gamma=1',n))
    xlabel('x')
    ylabel('y')
    count = count + 1;
end
figure(2)
sgtitle('Error Analysis for gamma = 1')
subplot(1,2,1)
plot(errMat(:,1), errMat(:,3))
title('n vs error')
xlabel('n')
ylabel('Error')
subplot(1,2,2)
plot(errMat(:,2), errMat(:,3))
title('l vs error')
xlabel('l')
ylabel('Error')

%%
%gamma = n + 1
figure(3)
count = 1;
errMat = [];%[n l error]
sgtitle('FEM Solutions')
for n = [16,32,64,128]
    l = 2;
    error = 1;
    [A,b] = popMatrices(n,n+1);  
    x0 = ones(n,1);
    M = eye(n);
    while error > 10^(-6)
        [sol,xs,ys,Vs,Hs] = gmres_matlab(A,b,l,x0, M, n);
        res = b - A*sol;
        error = norm(res)/n;
        l = l * 2;
    end
    errMat(end+1,:) = [n l error];
    
    x = linspace(0,1,n);
    subplot(2,2,count)
    plot(x,sol)
    title(sprintf('Solution with n = %d and gamma=%d',n,n+1))
    xlabel('x')
    ylabel('y')
    count = count + 1;
end
figure(4)
sgtitle('Error Analysis for gamma = n+1')
subplot(1,2,1)
plot(errMat(:,1), errMat(:,3))
title('n vs error')
xlabel('n')
ylabel('Error')
subplot(1,2,2)
plot(errMat(:,2), errMat(:,3))
title('l vs error')
xlabel('l')
ylabel('Error')
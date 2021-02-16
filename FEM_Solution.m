%gamma = 1
figure(1)
count = 1;
errMat = [];%[n l error]
sgtitle('FEM Solutions')
for n = [16,32,64,128]
    l = 2;
    error = 1;
    [A,b] = popMatrices(n,1);  
    x0 = zeros(n,1);
    M = eye(n);
    while error > 10^(-6)
        [sol,xs,ys,Vs,Hs] = gmres_matlab(A,b,l,x0, M, n);
        res = b - A*sol;
        error = norm(res)/n;
        l = l * 2;
        errMat(end+1,:) = [n l error];
    end    
    
    x = linspace(0,1,n);
    subplot(2,2,count)
    plot(x,sol)
    title(sprintf('Solution with n = %d and gamma=1',n))
    xlabel('x')
    ylabel('y')
    count = count + 1;
end
figure(2)
x = errMat(:,1);
y = errMat(:,2);
z = errMat(:,3);
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
trisurf(tri,x,y,z)
title('Error Analysis for gamma = 1')
xlabel('n')
ylabel('l')
zlabel('error')


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
    x0 = zeros(n,1);
    M = eye(n);
    while error > 10^(-6)
        [sol,xs,ys,Vs,Hs] = gmres_matlab(A,b,l,x0, M, n);
        res = b - A*sol;
        error = norm(res)/n;
        l = l * 2;
        errMat(end+1,:) = [n l error];
    end    
    
    x = linspace(0,1,n);
    subplot(2,2,count)
    plot(x,sol)
    title(sprintf('Solution with n = %d and gamma=%d',n,n+1))
    xlabel('x')
    ylabel('y')
    count = count + 1;
end
figure(4)
x = errMat(:,1);
y = errMat(:,2);
z = errMat(:,3);
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
trisurf(tri,x,y,z)
title('Error Analysis for gamma = n+1')
xlabel('n')
ylabel('l')
zlabel('error')
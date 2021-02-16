function[xs,ys,Vs,Hs]=gmres_matlab(A,b,l,x0,M,n)
r0=b - A*x0;
beta=norm(r0);
V(:,1)=r0/norm(r0);
H=[];%zeros(n,n);
for j=1:l
    W(:,j)=A*V(:,j);
    for i=1:j
        H(i,j)=dot(W(:,j),V(:,i));
        W(:,j)=H(i,j)*V(:,i);
    end
    H(j+1,j)=norm(W(:,j));
    if H(j+1,j)==0
        break;
    end
    V(:,j+1)=W(:,j)/H(j+1,j);
end
[rows,~]=size(H);
a=zeros(rows,1);
a(1)=beta;
ys=lsqlin(H,a);
for i=1:length(ys)
    xs(i)=V(i,1:length(ys))*ys;
end
Vs=V;
Hs=H;
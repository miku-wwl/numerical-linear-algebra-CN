function [k1,k2]=ex616( )
clc
format short e
A=[0.1 0.1 0.1 0.1 0.1 0.1; 0.0 0.2 0.1 0.1 0.1 0.1; 0.0 0.0 0.3 0.1 0.1 0.1; 
     0.0 0.0 0.0 0.4 0.1 0.1; 0.0 0.0 0.0 0.0 0.5 0.1; 0.0 0.0 0.0 0.0 0.0 0.6];
F=[3.2 1.5 1.0 0.6 0.3 0.1; 4.5 2.8 1.0 0.6 0.3 0.1; 6.0 4.3 2.5 0.6 0.3 0.1;
     7.7 6.0 4.2 2.3 0.3 0.1; 9.6 7.9 6.1 4.2 2.2 0.1; 11.7 10.0  8.2 6.3 4.3 2.2];
phiX=@(X) F-A*X-X*A';
EX=@(X) A*X+X*A'+X-F;
k1=0; X1=zeros(6);  
E=EX(X1);
res1(1)=norm(E,'fro');
tic
while (k1<=40)
    k1=k1+1;
    X=phiX(X1);
    E=EX(X);
    res=norm(E,'fro');
    if mod(k1,8)==0,
        res1(k1/8+1)=res;
    end
    X1=X;
end    
toc
%
k2=0; X1=zeros(6);  
E1=EX(X1);
res2(1)=norm(E1,'fro');
tic
while (k2<100)
    k2=k2+1;
    X2=phiX(X1);
    E2=EX(X2);
    a1=-trace((E1-E2)'*E2)/norm(E1-E2,'fro')^2;
    a2=1-a1;
    X=a1*X1+a2*X2;
    X1=X; E1=E2;
    res=norm(E1,'fro');
    if mod(k2,4)==0,
        res2(k2/4+1)=res;
    end
    if res<1.0e-7,
        r2=res; break;
    end
end    
toc
res1
res2
r2
X

function [k1,k2,r1,r2]=ex615( )
clc
format short e
A=[0.1 0.1 0.1 0.1 0.1 0.1; 0.0 0.2 0.1 0.1 0.1 0.1; 0.0 0.0 0.3 0.1 0.1 0.1; 
     0.0 0.0 0.0 0.4 0.1 0.1; 0.0 0.0 0.0 0.0 0.5 0.1; 0.0 0.0 0.0 0.0 0.0 0.6];
F=[2.2 1.5 1.0  0.6  0.3  0.1; 2.5 1.8 1.0 0.6 0.3 0.1;  3.0 2.3 1.5 0.6 0.3 0.1;
    3.7 3.0 2.2 1.3 0.3 0.1; 4.6 3.9 3.1 2.2 1.2 0.1; 5.7 5.0 4.2 3.3 2.3 1.2];
phiX=@(X) X-A*X-X*A'+F;
EX=@(X) A*X+X*A'-F;
k1=0; X1=zeros(6);  
E=EX(X1);
res1(1)=norm(E,'fro');
tic
while (k1<1000)
    k1=k1+1;
    X=phiX(X1);
    E=EX(X);
    res=norm(E,'fro');
    if mod(k1,10)==0,
        res1(k1/10+1)=res;
    end
    if res<1.0e-7,
        r1=res; break;
    end
    X1=X;
end    
toc
%
k2=0; X1=zeros(6);  
E1=EX(X1);
res2(1)=norm(E1,'fro');
tic
while (k2<1000)
    k2=k2+1;
    X2=phiX(X1);
    E2=EX(X2);
    a1=-trace((E1-E2)'*E2)/norm(E1-E2,'fro')^2;
    a2=1-a1;
    X=a1*X1+a2*X2;
    X1=X; E1=E2;
    res=norm(E1,'fro');
    if mod(k2,5)==0,
        res2(k2/5+1)=res;
    end
    if res<1.0e-7,
        r2=res; break;
    end
end    
toc
res1
res2
X

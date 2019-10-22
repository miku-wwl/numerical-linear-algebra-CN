function [k1,k2,t1,t2]=ex617(alpha)
clc
format short e
ni=3; m=100;  I=eye(ni*m); F=I;
A = cell(m,m);  B = cell(m,m);  Ii = -eye(ni);
AI = [4 -1 -1; -1  4 -1; -1 -1  4];  
BI = [4 -1 0; -1  4 -1; 0 -1  4]; 
for i=1:m
    for j=1:m
        A{i,j}=zeros(ni); B{i,j}=zeros(ni);
    end
end
for i=2:m-1
     A{i,i-1}=Ii; A{i,i}=AI; A{i,i+1}=Ii;
     B{i,i-1}=Ii; B{i,i}=BI; B{i,i+1}=Ii;
end
A{1,1}=AI;  A{1,2}=Ii; A{m,m-1}=Ii; A{m,m}=AI;
B{1,1}=BI;  B{1,2}=Ii; B{m,m-1}=Ii; B{m,m}=BI;
A=cell2mat(A); B=cell2mat(B);
C=I+2*inv(alpha*A-I);
D=I+2*inv(alpha*B-I);
Y0=-0.5*alpha*(C-I)*F*(D-I);
X0=Y0;  Iter=0;
Y0=-0.5*alpha*(C-I)*F*(D-I);
X1=Y0;  Iter=0;
phiX=@(X) C*X*D+Y0;
EX=@(X) A*X+X*B-F;
k1=0;  
E=EX(X1);
tic
while (k1<=1000)
    k1=k1+1;
    X=phiX(X1);
    E=EX(X);
    res=norm(E,'fro');
    if res<1.0e-10,
        break;
    end
    X1=X;
end    
t1=toc;
%
k2=0; X1=Y0;
E1=EX(X1);
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
    if res<1.0e-10,
       break;
    end
end    
t2=toc;



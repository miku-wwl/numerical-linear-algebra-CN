function [Iter,time]=ex611(alpha) 
%clc; 
ni=3; m=200;  I=eye(ni*m); F=ones(ni*m);
A = cell(m,m);  B = cell(m,m);  Ii = -eye(ni);
C = cell(m,m);  D = cell(m,m); 
AI = [4 -1 0;  -1  4 -1; 0 -1  4];   
BI =  [4 1 1;  1  4 1; 1 1  4]; 
CI=[4.5 -1 -1; -1 4.5 -1; -1 -1 4.5];
DI= [4 -1 -1;  -1  4 -1; -1 -1  4]; 
for i=1:m
    for j=1:m
        A{i,j}=zeros(ni); B{i,j}=zeros(ni);
        C{i,j}=zeros(ni); D{i,j}=zeros(ni);
    end
end
for i=2:m-1
     A{i,i-1}=Ii; A{i,i}=AI; A{i,i+1}=Ii;
     B{i,i-1}=Ii; B{i,i}=BI; B{i,i+1}=Ii;
     C{i,i-1}=Ii; C{i,i}=CI; C{i,i+1}=Ii;
     D{i,i-1}=Ii; D{i,i}=DI; D{i,i+1}=Ii;
end
A{1,1}=AI;  A{1,2}=Ii; A{m,m-1}=Ii; A{m,m}=AI;
B{1,1}=BI;  B{1,2}=Ii; B{m,m-1}=Ii; B{m,m}=BI;
C{1,1}=CI;  C{1,2}=Ii; C{m,m-1}=Ii; C{m,m}=CI;
D{1,1}=DI;  D{1,2}=Ii; D{m,m-1}=Ii; D{m,m}=DI;
A=cell2mat(A); B=cell2mat(B);
C=cell2mat(C); D=cell2mat(D);
A1=alpha*A+C; C1=alpha*A-C;
D1=alpha*D+B; B1=alpha*D-B;
rA1=inv(A1);   rD1=inv(D1);
P=rA1*C1;   Q=B1*rD1;
Y0=2*alpha*rA1*F*rD1';
X0=Y0; Iter=0;
tic
while(Iter<=500)
    X=P*X0*Q+Y0;
    Err=norm(X-X0,'fro');
    if(Err<=1.e-10)
       break;
    end
    X0=X;  Iter=Iter+1;
end
time=toc;
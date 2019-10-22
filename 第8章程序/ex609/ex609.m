function [k1,t1,k2,t2]=ex609(alpha) 
%clc;  
ni=3; m=100;  I=eye(ni*m); F=I;
A = cell(m,m);  B = cell(m,m);  Ii = -eye(ni);
AI = [12 1 1;  1  17 0; 1 0  15];  
BI = [9 1 1; 1 8 1; 1 1 7]; 
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
A=cell2mat(A);
B=cell2mat(B);
C=alpha*A+B; D=alpha*A-B;
rC=inv(C);   G=rC*D;
Y0=2*alpha*rC*F*rC';
S0=Y0; k1=0;
tic
while(k1<=500)
    S=S0+G^(2^k1)*S0*(G')^(2^k1);
    Err=norm(S-S0,'fro');
    if(Err<=1.e-10)
       break;
    end
    S0=S;  k1=k1+1;
end
k1=k1-1;
t1=toc;
X0=Y0; k2=0;
tic
while(k2<=500)
    X=G*X0*G'+Y0;
    Err=norm(X-X0,'fro');
    if(Err<=1.e-10)
       break;
    end
    X0=X;  k2=k2+1;
end
t2=toc;

%Àý 7.3
clear all
n=100; e=[1:n]'; e1=ones(n-1,1);
A=diag(e)-diag(e1,1)-diag(e1,-1);
[lam1,v1, k1]=mvpower(A,e,101);
[lam2,v2, k2]=mvpower(A,e,99);
[lam3,v3, k3]=mvpower(A,e,2);
[lam4,v4, k4]=mvpower(A,e,0);
d=eig(A);
Eigs=[d(1),d(2),d(end-1),d(end)],
Eigt=[lam4,lam3,lam2,lam1],
k=[k1,k2,k3,k4],


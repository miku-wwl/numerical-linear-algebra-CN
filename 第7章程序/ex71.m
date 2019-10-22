%Àý7.1
clear all
n=100; e=[1:n]'; e1=ones(n-1,1);
A=diag(e)-diag(e1,1)-diag(e1,-1);
[lam1,v, k]=mypower(A,e,1.e-6);
d=eig(A); 
d=max(abs(d))
lam1
k


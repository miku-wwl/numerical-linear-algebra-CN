%Àý7.2
clear all
n=100; e=[1:n]'; e1=ones(n-1,1);
A=diag(e)-diag(e1,1)-diag(e1,-1);
alpha=50;
[lam,v, k]=mopower(A,e,alpha);
d=eig(A); 
d=max(abs(d))
lam
k


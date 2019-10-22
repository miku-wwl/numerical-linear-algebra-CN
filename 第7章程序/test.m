clear all;
n=40000;
e=[2:n+1]';  e1=0.5*ones(n,1); e2=0.5*ones(n-1,1);
%A=spdiags([e1,e,-e1], -1:1, n, n);
%A1=diag(e)+diag(e2,-1)-diag(e2,1);
B=zeros(n,n);

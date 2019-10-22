%Àý 456
%clear all
n=10; e=6*ones(n,1); e1=ones(n-1,1); e2=ones(n-2,1); 
e3=ones(n-3,1); e4=ones(n-4,1);
A=diag(e)+diag(e1,1)-diag(e1,-1)-2*diag(e2,2)+2*diag(e2,-2);
%A=A+3*diag(e3,3)-3*diag(e3,-3)-2*diag(e4,4)+2*diag(e4,-4);
%A=[4,1,0,0;-1,4 ,1 ,0; 0, -1, 4, 1; 0, 0, -1,4];
v=rand(size(A,1),1);  %v=0.5*ones(size(A,1),1); 
[mu1,u1,Vm]=jacobidavidson(A,v,1.e-12);
mu1,
% [mu2,u2,Vm]=jacobidavi(A,v,8);
% mu2,
d=eig(A); 
[s,t]=max(d);
s=d(t)







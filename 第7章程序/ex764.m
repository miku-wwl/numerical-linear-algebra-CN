%Àý7.6.4-ex764.m
clear all; tic;
n=1000; A=zeros(n,n);
for i=1:n,
    for j=1:i
        A(i,j)=j;  A(j,i)=A(i,j);
    end
end
v=rand(size(A,1),1);  
%alpha=18; B=inv(A-alpha*eye(n));
[mu,u,Vm]=jacobidavidson(A,v);  %Vm'*Vm
%mu=alpha+1/mu
mu
err1=norm(A*u-mu*u)
d=eig(A);
[s,t]=max(abs(d));
err2=norm(mu-d(t))
toc







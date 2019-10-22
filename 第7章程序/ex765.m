%Àý7.6.5-ex765.m
clear all; tic;
n=1000; A=zeros(n,n);
for i=1:n,
    for j=1:i
        A(i,j)=j;  A(j,i)=A(i,j);
    end
end
v=rand(n,1); I=eye(n);
alpha=18; B=(A-alpha*I)\I;
[mu,u,Vm]=jacobidavidson(B,v);  %Vm'*Vm
lambda=alpha+1/mu
err1=norm(A*u-lambda*u)
d=eig(A);
err2=norm(lambda-d(925))
toc



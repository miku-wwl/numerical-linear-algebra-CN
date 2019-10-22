%Àý7.6.2-ex762.m
clear all
n=2000; m=5; A=zeros(n,n);
for i=1:n 
    for j=1:i
        A(i,j)=j;  A(j,i)=A(i,j);
    end
end
ev=eig(A); ev=sort(ev,'descend'); lam=ev(1:m);
[M,V,iter] =irest_Lanczos(A,m); iter
mu=diag(M(1:m,1:m))
err1=norm(mu-lam)
err2=norm(A*V(:,1)-M(1,1)*V(:,1))






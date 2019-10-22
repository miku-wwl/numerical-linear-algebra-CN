function [A]=hessenb(A)
%��������Householder�任��n�׾���AΪ��Hessenberg����.
%���ú���: mhouse.m
n=size(A,1);  %Q=eye(n);
for k=1:(n-2)
    x=A(k+1:n,k); [v,beta]=mhouse(x);
    H=(eye(length(v))-beta*v*v');
    A(k+1:n,1:n)=H*A(k+1:n,1:n);
    A(1:n,k+1:n)=A(1:n,k+1:n)*H;
    %Q=Q*blkdiag(eye(k),H);
end
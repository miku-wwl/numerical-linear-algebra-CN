%Cholesky分解法程序--mschol.m
function [x,L,D]=mschol(A,b)
%用Cholesky分解法解对称正定方程组 Ax=b
%输入矩阵A, 右端项b; 输出解向量x,单位下三角阵L,对角阵D
n=size(A,1); D=zeros(1,n); L=eye(n,n);
U(1,:)=A(1,:); L(2:n,1)=U(1,2:n)/U(1,1); %Cholesky分解
for k=2:n
   U(k,k:n)=A(k,k:n)-L(k,1:k-1)*U(1:k-1,k:n);
   L(k+1:n,k)=U(k,k+1:n)/U(k,k);
end
%求解下三角方程组Ly=b(向前消去法)
y=zeros(n,1);  y(1)=b(1);
for k=2:n,
    y(k)=b(k)-L(k,1:k-1)*y([1:k-1]);
end
%求解对角方程组Dz=y
D=diag(diag(U));
for k=1:n, 
    z(k)=y(k)/D(k,k); 
end
%求解上三角方程组L'x=z(回代法)
x=zeros(n,1);   U=L';  x(n)=z(n);
for k=(n-1):-1:1,
    x(k)=z(k)-U(k,k+1:n)*x(k+1:n);
end
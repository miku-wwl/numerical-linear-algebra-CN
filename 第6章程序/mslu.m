%顺序LU分解法程序--mslu.m
function [x,A]=mslu(A,b)
%顺序LU分解A=LU, A为系数矩阵, b为右端向量.
%x为解向量, L和U分别存放在A的严格下三角和上三角部分
n=length(b);
for k=1:n  %顺序LU分解
   A(k:n,k)=A(k:n,k)-A(k:n,1:k-1)*A(1:k-1,k);
   A(k+1:n,k)=A(k+1:n,k)/A(k,k); %乘子向量
   A(k,k+1:n)=A(k,k+1:n)-A(k,1:k-1)*A(1:k-1,k+1:n);
end
y=zeros(n,1);
for k=1:n, %解下三角矩阵Ly=b
   y(k)=b(k)-A(k,1:k-1)*y(1:k-1);
end
x=zeros(n,1);
for k=n:-1:1, %解上三角方程组Ux=y
   x(k)=(y(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
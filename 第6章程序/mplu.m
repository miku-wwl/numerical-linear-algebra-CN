%列主元LU分解法程序--mplu.m
function [x,A,P]=mplu(A,b)
%列主元LU分解PA=LU, A为系数矩阵, b为右端向量,
%x返回解向量, L和U分别存放在A的严格下三角和
%上三角部分,P返回选主元时记录行交换的置换阵
n=length(b);
P=eye(n); %P记录选择主元时候所进行的行变换
for k=1:n  %列主元LU分解
   A(k:n,k)=A(k:n,k)-A(k:n,1:k-1)*A(1:k-1,k);
   [s,m]=max(abs(A(k:n,k)));   %选列主元
   m=m+k-1;
   if m~=k
      A([k m],:)=A([m k],:);
      P([k m],:)=P([m k],:);
      %b([k m],:)=b([m k],:);
   end
   A(k+1:n,k)=A(k+1:n,k)/A(k,k);
   A(k,k+1:n)=A(k,k+1:n)-A(k,1:k-1)*A(1:k-1,k+1:n);
end
b=P*b; y=zeros(n,1);
for k=1:n, %解下三角矩阵Ly=b
   y(k)=b(k)-A(k,1:k-1)*y(1:k-1);
end
x=zeros(n,1);
for k=n:-1:1, %解上三角方程组Ux=y
   x(k)=(y(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
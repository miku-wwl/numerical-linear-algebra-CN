%顺序Gauss消去法程序-msgauss.m
function x=msgauss(A,b,flag)
%输入: A为系数矩阵, b为右端项, 若flag=0,则不显示中
%间过程, 否则显示中间过程, 默认为0. 输出: x为解向量
if nargin<3, flag=0; end
n=length(b);
for k=1:(n-1)  %消元过程
   m=A(k+1:n,k)/A(k,k);
   A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-m*A(k,k+1:n);
   b(k+1:n)=b(k+1:n)-m*b(k);
   A(k+1:n,k)=zeros(n-k,1);
   if flag~=0, Ab=[A,b], end
end
x=zeros(n,1); %回代过程
x(n)=b(n)/A(n,n);
for k=n-1:-1:1
   x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
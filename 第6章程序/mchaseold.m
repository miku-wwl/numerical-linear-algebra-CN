%追赶法程序--mchase.m
function [x]=mchaseold(a,b,c,f)
%用追赶法解三对角方程组 Ax=f
%输入: a=[0, a(2),a(3),...,a(n)] 为A的下对角线，
%b=[b(1),b(2),...,b(n-1),b(n)] 为A的主对角线，
%c=[c(1),c(2),...,c(n-1),0] 为A的次上对角线, 
%f为右端向量. 输出: 解向量x
n=length(b);
l(1)=b(1); y(1)=f(1)/l(1);
for k=2:n  
    u(k-1)=c(k-1)/l(k-1);
    l(k)=b(k)-a(k)*u(k-1);
    y(k)=(f(k)-a(k)*y(k-1))/l(k);
end    
x(n)=y(n);
for k=n-1:-1:1
    x(k)=y(k)-u(k)*x(k+1);
end
    
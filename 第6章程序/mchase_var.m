%变参数追赶法程序--mchase_var.m
function [x]=mchase_var(a,b,c,f,l1,u1)
%用变参数追赶法解三对角方程组 Ax=f
%输入: a为A的下对角线,b为A的主对角%线 c为A的
%次上对角线, f为右端向量. l1,u1满足l1*u1+1不为0
%输出: 解向量x
n=length(b);  x=zeros(n,1); 
u=zeros(n,1); l=zeros(n,1); 
l(1)=l1; u(1)=u1; 
d(1)=b(1)/(l(1)*u(1)+1); g(1)=f(1)/d(1);
for k=2:n  
    u(k)=c(k-1)/d(k-1);
    d(k)=b(k)-a(k)*u(k);
    l(k)=a(k)/d(k);  g(k)=f(k)/d(k);
end    
s(n)=1+u(n)*l(n); t(n)=g(n);
for k=n-1:-1:1
    s(k)=1+u(k)*s(k+1)*l(k); 
    t(k)=s(k+1)*g(k)-u(k+1)*t(k+1);
end
x(1)=t(1)/s(1); y(1)=u(1)*x(1);
for k=1:n
    y(k+1)=g(k)-l(k)*y(k);
end
x(n)=y(n+1);
for k=n-1:-1:2
    x(k)=y(k+1)-u(k+1)*x(k+1);
end

    
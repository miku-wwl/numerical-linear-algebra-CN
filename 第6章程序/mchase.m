%追赶法程序--mchase.m
function [f]=mchase(a,b,c,f)
%用追赶法解三对角方程组 Ax=f
%输入: a为A的下对角线,b为A的主对角
%线 c为A的次上对角线, f为右端向量. 
%输出: 解向量f (LU分解中的l(k),u(k)存放在b(k),
%c(k)的位置, y(k)和x(k)先后存放在d(k)的位置)
n=length(b); f(1)=f(1)/b(1);
for k=2:n  
    c(k-1)=c(k-1)/b(k-1);
    b(k)=b(k)-a(k)*c(k-1);
    f(k)=(f(k)-a(k)*f(k-1))/b(k);
end    
for k=n-1:-1:1
    f(k)=f(k)-c(k)*f(k+1);
end
    
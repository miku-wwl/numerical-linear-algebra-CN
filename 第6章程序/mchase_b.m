%块追赶法程序--mchase_b.m
function [x]=mchase_b(Ai,Bi,Ci,fi,m)
%用追赶法解三对角方程组 Ax=f
%输入: Ai为A的次下对角块,Bi为A的主对角块,
%Ci为A的次上对角块, fi为右端向量. 
%输出: 解向量x (LU分解中的l(k),u(k)存放在b(k),
%c(k)的位置, y(k)和x(k)先后存放在d(k)的位置)
x=cell(m,1);
L{1}=Bi{1}; y{1}=L{1}\fi{1};
for k=2:m  
    U{k-1}=L{k-1}\Ci{k-1};
    L{k}=Bi{k}-Ai{k}*U{k-1};
    y{k}=L{k}\(fi{k}-Ai{k}*y{k-1});
end  
x{m}=y{m};
for k=m-1:-1:1
    x{k}=y{k}-U{k}*x{k+1};
end
    
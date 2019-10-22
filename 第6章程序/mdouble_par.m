%双参数法程序--mdouble_par.m
function [x]=mdouble_par(Ai,Bi,Ci,fi,m)
%用双参数法解三对角方程组 Ax=f
%输入: Ai为A的次下对角块,Bi为A的主对角块,
%Ci为A的次上对角块, fi为右端向量.
%输出: 解向量x
x=cell(m,1); s=cell(m,1); T=cell(m,1);
s{1}=zeros(3,1); s{2}=Ci{1}\fi{1}; 
T{1}=eye(3); T{2}=-Ci{1}\Bi{1};
for k=2:m-1 
    s{k+1}=Ci{k}\(fi{k}-Ai{k}*s{k-1}-Bi{k}*s{k});
    T{k+1}=-Ci{k}\(Ai{k}*T{k-1}+Bi{k}*T{k});
end  
x{1}=(Ai{m}*T{m-1}+Bi{m}*T{m})\(fi{m}-Ai{m}*s{m-1}-Bi{m}*s{m});
for k=2:m
    x{k}=s{k}+T{k}*x{1};
end
    
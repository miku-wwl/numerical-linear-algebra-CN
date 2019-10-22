%块追赶法程序--mchase_block.m
function [fi]=mchase_block(Ai,Bi,Ci,fi,m)
%用块追赶法解块三对角方程组 Ax=f
%输入: Ai为A的次下对角块,Bi为A的主对角块,
%Ci为A的次上对角块, fi为右端向量,m为分块数. 
%输出: 解向量fi,其中LU分解的L{k},U{k}存放在Bi{k},
%Ci{k}的位置,y{k}和x{k}先后存放在fi{k}的位置
fi{1}=Bi{1}\fi{1};
for k=2:m  
    Ci{k-1}=Bi{k-1}\Ci{k-1};
    Bi{k}=Bi{k}-Ai{k}*Ci{k-1};
    fi{k}=Bi{k}\(fi{k}-Ai{k}*fi{k-1});
end  
for k=m-1:-1:1
    fi{k}=fi{k}-Ci{k}*fi{k+1};
end
    
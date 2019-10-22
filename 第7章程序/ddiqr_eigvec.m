function [Lam,V,iter,ki]=ddiqr_eigvec(A,tol)
%用双重位移隐式QR方法求实方阵的全部特征值和相应的特征向量.
%输入: n阶实方阵A, 控制精度tol(默认是1.e-5)
%输出: 迭代次数Iter, A的全部特征值Lam和特征向量V
if nargin<2, tol=1e-5; end
n=size(A,1);  x=rand(n,1); %x=ones(n,1);
Lam=zeros(n,1); V=zeros(n);  
[A,Q]=mhessen(A); %调用上Hessenberg化程序
[iter,lambda]=ddiqr_eig(A,tol); %调用双重位移隐式QR方法求全部特征值
for i=1:n
    [lam,v,k]=mvpower(A,x,lambda(i)); %调用反幂法程序
    V(:,i)=v; ki(i)=k; Lam(i)=lam;
end
V=Q*V;    %V的每一列为特征向量




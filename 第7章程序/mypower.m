%幂法程序--mypower.m
function [lam,v,k]=mypower(A,x,tol,N)
%本程序用幂法计算矩阵的模最大特征值和对应的特征
%向量,其中输入参数A为n阶方阵,x为初始向量, tol控制
%精度,N为最大迭代次数. 输出参数lam为按模最大 
%的特征值,v为对应的特征向量, k为迭代次数.
if nargin<4, N=1000; end
if nargin<3, tol=1e-6;end
m=0; k=0; err=1;
while(k<N)
   v=A*x;
   [m1,t]=max(abs(v));
   m1=v(t);   x=v/m1;
   err=abs(m1-m);
   if err<tol, break; end
   m=m1;  k=k+1;
end
lam=m1; v=x;
%反幂法程序-mvpower.m
function [lam,v,k]=mvpower(A,x,alpha,tol,N)
%本程序用反幂法计算矩阵与alpha最接近的特征值和对应的特征向量,
%输入:A为n阶方阵,x为初始向量,tol为精度,N为最大迭代数,alpha为某个常数.
%输出:lam返回与alpha最接近的特征值,v返回对应的特征向量,k返回迭代次数
if nargin<5, N=500; end
if nargin<4,tol=1e-5;end
m=0.5; k=0; err=1;
A=A-alpha*eye(length(x));
[L,U,P]=lu(A);
while (k<N)
   [m1,t]=max(abs(x));
   m1=x(t); v=x/m1;
   z=L\(P*v);  x=U\z;
   err=abs(1/m1-1/m);
   if err<=tol,break;end
   k=k+1;  m=m1;
end
lam=alpha+1/m;
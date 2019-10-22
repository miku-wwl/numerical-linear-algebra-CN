%原点位移幂法程序-mopowerp.m
function [lam,v,k]=mopower(A,x,alpha,tol,N)
%本程序用原点位移幂法求矩阵的模最大特征值和对应的特征向量,
%输入参数A为n阶方阵, x为初始向量, tol为控制精度, N最大迭代
%次数, alpha为原点位移参数. 输出参数lam返回按模最大的特征
%值, v返回对应的特征向量, k返回迭代次数.
if nargin<5,N=1000;end
if nargin<4,tol=1e-6;end
m=0; k=0; err=1;
A=A-alpha*eye(length(x));
while(k<N)
    v=A*x;
    [m1,t]=max(abs(v));
    m1=v(t);   x=v/m1;
    err=abs(m1-m); 
    if err<tol, break; end
    m=m1;  k=k+1;
end
lam=m1+alpha; v=x;
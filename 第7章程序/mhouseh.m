%程序--mhouseh.m
function [H,a]=mhouseh(x)
%用途: 对于向量x,构造Householder变换矩阵H,使得Hx=(-a,0,...,0)'
%格式: function H=mhouseh(x)
%x为输入列向量, H返回Householder变换矩阵
n=length(x); I=eye(n); sn=sign(x(1));
if sn==0, sn=1; end
z=x(2:n);
if(norm(z,inf)==0), H=I; return; end
a=sn*norm(x);
u=x; u(1)=u(1)+a;
rho=a*(a+x(1));
H=I-1.0/rho*u*u';
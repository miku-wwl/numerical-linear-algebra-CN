%����--mhouseh.m
function [H,a]=mhouseh(x)
%��;: ��������x,����Householder�任����H,ʹ��Hx=(-a,0,...,0)'
%��ʽ: function H=mhouseh(x)
%xΪ����������, H����Householder�任����
n=length(x); I=eye(n); sn=sign(x(1));
if sn==0, sn=1; end
z=x(2:n);
if(norm(z,inf)==0), H=I; return; end
a=sn*norm(x);
u=x; u(1)=u(1)+a;
rho=a*(a+x(1));
H=I-1.0/rho*u*u';
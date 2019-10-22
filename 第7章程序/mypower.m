%�ݷ�����--mypower.m
function [lam,v,k]=mypower(A,x,tol,N)
%���������ݷ���������ģ�������ֵ�Ͷ�Ӧ������
%����,�����������AΪn�׷���,xΪ��ʼ����, tol����
%����,NΪ����������. �������lamΪ��ģ��� 
%������ֵ,vΪ��Ӧ����������, kΪ��������.
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
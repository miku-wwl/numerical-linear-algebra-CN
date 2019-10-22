%���ݷ�����-mvpower.m
function [lam,v,k]=mvpower(A,x,alpha,tol,N)
%�������÷��ݷ����������alpha��ӽ�������ֵ�Ͷ�Ӧ����������,
%����:AΪn�׷���,xΪ��ʼ����,tolΪ����,NΪ��������,alphaΪĳ������.
%���:lam������alpha��ӽ�������ֵ,v���ض�Ӧ����������,k���ص�������
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
%ԭ��λ���ݷ�����-mopowerp.m
function [lam,v,k]=mopower(A,x,alpha,tol,N)
%��������ԭ��λ���ݷ�������ģ�������ֵ�Ͷ�Ӧ����������,
%�������AΪn�׷���, xΪ��ʼ����, tolΪ���ƾ���, N������
%����, alphaΪԭ��λ�Ʋ���. �������lam���ذ�ģ��������
%ֵ, v���ض�Ӧ����������, k���ص�������.
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
%Givens-Householder��������-givens_househ.m
function [lambda,k]=givens_househ(T,m,tol,max_it)
%��������Givens-Householder��ʵ�Գ����ԽǾ���ĵ�m������ֵ.
%�������:Ϊn��ʵ�ԳƷ���, m����ֵ���,tolΪ�������.
%�������:lambda��m������ֵ,kΪ���㾫�ȵ�������
if nargin<4, max_it=500; end
if nargin<3, tol=1e-6; end
n=size(T,1);  k=0;
b=norm(T,inf); a=-b;
alpha=diag(T); beta=diag(T,1);
while(k<max_it)
    k=k+1; c=(a+b)/2; 
    q=qkfun(c,alpha,beta);
    tt=find(q>=0); 
    sn=length(tt); 
    if(sn>=m), a=c; else, b=c; end
    if abs(b-a)<=tol,  
        lambda=c; break; 
    end
end
   
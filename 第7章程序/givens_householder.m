%Givens-Householder��������-givens_householder.m
function [Lam,Q]=givens_householder(A,tol)
%��������Givens-Householder��ʵ�Գƾ���A�ĵ�m������ֵ.
%�������:AΪn��ʵ�ԳƷ���, m����ֵ���,tolΪ�������.
%�������:lambda��m������ֵ,kΪ���㾫�ȵ�������
if nargin<4, max_it=500; end
if nargin<3, tol=1e-6; end
n=size(A,1);  k=0; Lam=zeros(n,1);
[T,Q]=mhessen(A); b=norm(T,inf); a=-b;
alpha=diag(T); beta=diag(T,1);
for i=1:n
    while(1)
        c=(a+b)/2;
        q=qkfun(c,alpha,beta);
        tt=find(q>=0); 
        sn=length(tt); 
        if(sn>=i), a=c; else, b=c; end
        if abs(b-a)<=tol,
            Lam(i)=c; break; 
        end
    end
end
   
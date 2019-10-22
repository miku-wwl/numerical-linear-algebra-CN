function [Lam,V,iter,ki]=ddiqr_eigvec(A,tol)
%��˫��λ����ʽQR������ʵ�����ȫ������ֵ����Ӧ����������.
%����: n��ʵ����A, ���ƾ���tol(Ĭ����1.e-5)
%���: ��������Iter, A��ȫ������ֵLam����������V
if nargin<2, tol=1e-5; end
n=size(A,1);  x=rand(n,1); %x=ones(n,1);
Lam=zeros(n,1); V=zeros(n);  
[A,Q]=mhessen(A); %������Hessenberg������
[iter,lambda]=ddiqr_eig(A,tol); %����˫��λ����ʽQR������ȫ������ֵ
for i=1:n
    [lam,v,k]=mvpower(A,x,lambda(i)); %���÷��ݷ�����
    V(:,i)=v; ki(i)=k; Lam(i)=lam;
end
V=Q*V;    %V��ÿһ��Ϊ��������




function [mu,U,lambda] = Class_Lanczos(A,v,k)
%����Lanczos����.����Գƾ���A,��ʼ����v���ӿռ�ά��k
%���: lambda��U�ֱ���A��k������ֵ����Ӧ����������.
tol=1.e-12;
v=v/norm(v); V(:,1)=v;
for i=1:k
    [V,T,beta]=Lanczos2(A,v,i);
    [Y,Mu]=eig(T);
    mu{i}=diag(Mu);
    U=V*Y;
    if beta*abs(Y(i,i))<tol
       k=i; break; 
    end
end
lambda=mu{k};




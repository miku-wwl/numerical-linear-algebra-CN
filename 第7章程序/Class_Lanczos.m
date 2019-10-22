function [mu,U,lambda] = Class_Lanczos(A,v,k)
%经典Lanczos方法.输入对称矩阵A,初始向量v和子空间维数k
%输出: lambda和U分别是A的k个特征值及相应的特征向量.
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




function [kappa]=cond_inf(A)
%输入: 矩阵A; 输出: 矩阵A无穷范数条件数的估计值
n=size(A,1);
[L,U] = lu(A); %列主元LU分解PA=LU
%求解U'x=b
beta(1)=0; x=zeros(n,1);
b(1)=1; x(1)=b(1)/U(1,1);
for k=2:n
    beta(k)=U(1:k-1,k)'*x(1:k-1); 
    a=beta(k); b(k)=-sign(a);
    x(k)=(b(k)-beta(k))/U(k,k);
end
dw=L'\x; z=L\dw; y=U\z;
h2=norm(y,'inf')/norm(dw,'inf');
h1=norm(A,'inf'); 
kappa=h1*h2; %矩阵A的无穷范数条件数












        
    
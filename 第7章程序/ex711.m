%例7.11
n=15; tol=1.e-12; max_it=1000;
A=zeros(n); lambda=zeros(n,1); iter=zeros(n,1); V=zeros(n,n);
for i=1:n,
    for j=1:i
        A(i,j)=j;  A(j,i)=A(i,j);
    end
end
[Lam,V,ki]=ghvector(A,tol,max_it);
D=[sort(Lam),   eig(A)],  %显示计算结果
%V   %V的每一列为特征向量
norm(A*V-V*diag(Lam),inf)  %验证计算是否正确


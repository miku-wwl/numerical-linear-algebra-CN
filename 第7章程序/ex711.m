%��7.11
n=15; tol=1.e-12; max_it=1000;
A=zeros(n); lambda=zeros(n,1); iter=zeros(n,1); V=zeros(n,n);
for i=1:n,
    for j=1:i
        A(i,j)=j;  A(j,i)=A(i,j);
    end
end
[Lam,V,ki]=ghvector(A,tol,max_it);
D=[sort(Lam),   eig(A)],  %��ʾ������
%V   %V��ÿһ��Ϊ��������
norm(A*V-V*diag(Lam),inf)  %��֤�����Ƿ���ȷ


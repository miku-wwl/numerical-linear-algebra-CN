%Cholesky�ֽⷨ����--mschol.m
function [x,L,D]=mschol(A,b)
%��Cholesky�ֽⷨ��Գ����������� Ax=b
%�������A, �Ҷ���b; ���������x,��λ��������L,�Խ���D
n=size(A,1); D=zeros(1,n); L=eye(n,n);
U(1,:)=A(1,:); L(2:n,1)=U(1,2:n)/U(1,1); %Cholesky�ֽ�
for k=2:n
   U(k,k:n)=A(k,k:n)-L(k,1:k-1)*U(1:k-1,k:n);
   L(k+1:n,k)=U(k,k+1:n)/U(k,k);
end
%��������Ƿ�����Ly=b(��ǰ��ȥ��)
y=zeros(n,1);  y(1)=b(1);
for k=2:n,
    y(k)=b(k)-L(k,1:k-1)*y([1:k-1]);
end
%���ԽǷ�����Dz=y
D=diag(diag(U));
for k=1:n, 
    z(k)=y(k)/D(k,k); 
end
%��������Ƿ�����L'x=z(�ش���)
x=zeros(n,1);   U=L';  x(n)=z(n);
for k=(n-1):-1:1,
    x(k)=z(k)-U(k,k+1:n)*x(k+1:n);
end
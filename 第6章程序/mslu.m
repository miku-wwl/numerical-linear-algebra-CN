%˳��LU�ֽⷨ����--mslu.m
function [x,A]=mslu(A,b)
%˳��LU�ֽ�A=LU, AΪϵ������, bΪ�Ҷ�����.
%xΪ������, L��U�ֱ�����A���ϸ������Ǻ������ǲ���
n=length(b);
for k=1:n  %˳��LU�ֽ�
   A(k:n,k)=A(k:n,k)-A(k:n,1:k-1)*A(1:k-1,k);
   A(k+1:n,k)=A(k+1:n,k)/A(k,k); %��������
   A(k,k+1:n)=A(k,k+1:n)-A(k,1:k-1)*A(1:k-1,k+1:n);
end
y=zeros(n,1);
for k=1:n, %�������Ǿ���Ly=b
   y(k)=b(k)-A(k,1:k-1)*y(1:k-1);
end
x=zeros(n,1);
for k=n:-1:1, %�������Ƿ�����Ux=y
   x(k)=(y(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
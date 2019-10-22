%����ԪLU�ֽⷨ����--mplu.m
function [x,A,P]=mplu(A,b)
%����ԪLU�ֽ�PA=LU, AΪϵ������, bΪ�Ҷ�����,
%x���ؽ�����, L��U�ֱ�����A���ϸ������Ǻ�
%�����ǲ���,P����ѡ��Ԫʱ��¼�н������û���
n=length(b);
P=eye(n); %P��¼ѡ����Ԫʱ�������е��б任
for k=1:n  %����ԪLU�ֽ�
   A(k:n,k)=A(k:n,k)-A(k:n,1:k-1)*A(1:k-1,k);
   [s,m]=max(abs(A(k:n,k)));   %ѡ����Ԫ
   m=m+k-1;
   if m~=k
      A([k m],:)=A([m k],:);
      P([k m],:)=P([m k],:);
      %b([k m],:)=b([m k],:);
   end
   A(k+1:n,k)=A(k+1:n,k)/A(k,k);
   A(k,k+1:n)=A(k,k+1:n)-A(k,1:k-1)*A(1:k-1,k+1:n);
end
b=P*b; y=zeros(n,1);
for k=1:n, %�������Ǿ���Ly=b
   y(k)=b(k)-A(k,1:k-1)*y(1:k-1);
end
x=zeros(n,1);
for k=n:-1:1, %�������Ƿ�����Ux=y
   x(k)=(y(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
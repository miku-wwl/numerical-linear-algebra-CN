function [Q,A]=mhessen_qrtran(A)
%����: ��Givens�任����Hessenberg����A����QR�ֽ�
%����: n����Hessenberg����A, ����A(i+1,i)=0, i>2.
%���: ��������Q�������Ǿ���A: QA=A.
n=size(A,1); Q=eye(n);
for i=1:n-1
   xi=A(i,i);  xk=A(i+1,i);
   if xk~=0
      d=sqrt(xi^2+xk^2);
      c=xi/d;   s=xk/d;
      G=[c, s; -s, c];
      A(i:i+1,i:n)=G*A(i:i+1,i:n);
      Q(1:n,i:i+1)=Q(1:n,i:i+1)*G';
   end
end

function [A]=hessen_qrtran(A,m)
%����������n����Hessenberg����A,��Givens�任�������Ͻ�m�����ӿ�
%����QR�ֽ�,�������Ʊ任,�������任�����Hessenberg����A.
Q=eye(m);
for i=1:m-1
   xi=A(i,i);  xk=A(i+1,i);
   if xk~=0
      d=sqrt(xi^2+xk^2);
      c=xi/d;   s=xk/d;
      G=[c, s; -s, c];
      A(i:i+1,i:m)=G*A(i:i+1,i:m);
      Q(1:m,i:i+1)=Q(1:m,i:i+1)*G';
   end
end
%Q*A,  %��֤Q*R=A
A(1:m,1:m)=A(1:m,1:m)*Q;

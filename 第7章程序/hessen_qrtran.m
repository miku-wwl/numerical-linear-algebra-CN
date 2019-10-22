function [A]=hessen_qrtran(A,m)
%本程序输入n阶上Hessenberg矩阵A,用Givens变换对其左上角m阶主子块
%进行QR分解,再做相似变换,最后输出变换后的上Hessenberg矩阵A.
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
%Q*A,  %验证Q*R=A
A(1:m,1:m)=A(1:m,1:m)*Q;

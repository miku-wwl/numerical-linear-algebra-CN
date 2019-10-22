function [Q,A]=mhessen_qrtran(A)
%功能: 用Givens变换对上Hessenberg矩阵A进行QR分解
%输入: n阶上Hessenberg矩阵A, 其中A(i+1,i)=0, i>2.
%输出: 正交矩阵Q和上三角矩阵A: QA=A.
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

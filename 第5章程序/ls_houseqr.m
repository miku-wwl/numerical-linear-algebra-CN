function [x]=ls_houseqr(A,b)
%�㷨5.2, Householder QR�ֽ���min||Ax-b||
[n]=size(A,2);
[A]=house_qr([A,b]);
x=triu(A(1:n,1:n))\A(1:n,n+1);


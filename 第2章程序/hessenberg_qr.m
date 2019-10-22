function [A]=hessenberg_qr(A)
%上Hessenberg矩阵的QR分解
[n]=size(A, 2);
for k=1:n-1
    [c,s]=givens(A(k,k), A(k+1,k));
    A(k:k+1, k:n)=[c, s; -s, c]*A(k:k+1, k:n);
end
    
function [A]=givens_qr(A)
%Givens变换QR分解,不显式计算并存储正交阵Q
[m,n]=size(A);
for k=1:n
    for i=m:-1:k+1
        [c,s]=r_givens(A(i-1,k), A(i,k));
        A(i-1:i, k:n)=[c, s; -s, c]*A(i-1:i, k:n);
    end
end
    
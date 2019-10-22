function [Q,R]=Qhouse_qr(A)
%Householder QR分解, 显式计算并存储Q
[m,n]=size(A); Q=eye(m);
for k=1:n
    if k<m
        [v,beta]=r_house(A(k:m,k)); 
        H=eye(m-k+1)-beta*v*v';
        A(k:m,k:n)=H*A(k:m,k:n); 
        Q=Q*blkdiag(eye(k-1),H);   
    end
end
R=triu(A);

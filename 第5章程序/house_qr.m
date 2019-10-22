function [A,d]=house_qr(A)
%À„∑®2.3, Householder QR∑÷Ω‚
[m,n]=size(A); %Q=eye(n);
for k=1:n
    if k<m
        [v,beta]=r_house(A(k:m,k));
        A(k:m,k:n)=(eye(m-k+1)-beta*v*v')*A(k:m,k:n); 
        d(k)=beta;  A(k+1:m,k)=v(2:m-k+1);     
    end
end
%R=triu(A);

function [Q,R]=Qgivens_qr(A)
%Givens�任QR�ֽ�,��ʽ���㲢�洢������Q
[m,n]=size(A); Q=eye(m);
for k=1:n
    for i=m:-1:k+1
        [c,s]=r_givens(A(i-1,k), A(i,k));
        A(i-1:i, k:n)=[c, s; -s, c]*A(i-1:i, k:n);
        Q(:,i-1:i)=Q(:,i-1:i)*[c, s; -s, c]';
    end
end
R=triu(A);
    
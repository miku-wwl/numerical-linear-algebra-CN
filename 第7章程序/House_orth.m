%Householder 正交化程序-House_orth.m
function [Q,R]=House_orth(X)
[n,m]=size(X); %矩阵的行和列
R(:,1)=X(:,1);  E=eye(n); Y=cell(m,1);
for k=1:m 
    if k>1
        z=X(:,k); 
        for i=1:k-1
            v=Y{i};  %取出来
            H=blkdiag(eye(i-1),eye(n-i+1)-b(i)*v*v');
            z=H*z;
        end
        R(:,k)=z;
    end
    [v,beta]=r_house(R(k:n,k)); 
    H=blkdiag(eye(k-1),eye(n-k+1)-beta*v*v');
    R(:,k)=H*R(:,k);
    b(k)=beta;  Y{k}=v;  %存起来
    z=E(:,k);
    for i=k:-1:1
        v=Y{i};  %取出来
        H=blkdiag(eye(i-1),eye(n-i+1)-b(i)*v*v');
        z=H*z;
    end
    Q(:,k)=z;
end
R=R(1:m,:);

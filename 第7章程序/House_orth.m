%Householder ����������-House_orth.m
function [Q,R]=House_orth(X)
[n,m]=size(X); %������к���
R(:,1)=X(:,1);  E=eye(n); Y=cell(m,1);
for k=1:m 
    if k>1
        z=X(:,k); 
        for i=1:k-1
            v=Y{i};  %ȡ����
            H=blkdiag(eye(i-1),eye(n-i+1)-b(i)*v*v');
            z=H*z;
        end
        R(:,k)=z;
    end
    [v,beta]=r_house(R(k:n,k)); 
    H=blkdiag(eye(k-1),eye(n-k+1)-beta*v*v');
    R(:,k)=H*R(:,k);
    b(k)=beta;  Y{k}=v;  %������
    z=E(:,k);
    for i=k:-1:1
        v=Y{i};  %ȡ����
        H=blkdiag(eye(i-1),eye(n-i+1)-b(i)*v*v');
        z=H*z;
    end
    Q(:,k)=z;
end
R=R(1:m,:);

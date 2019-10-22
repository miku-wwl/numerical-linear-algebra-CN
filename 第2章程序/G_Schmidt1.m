%��׼Gram-Schmidt ����������-G_Schmidt1.m
function [Q,R]=G_Schmidt1(X)
[m]=size(X,2); %�������
R=zeros(m);  R(1,1)=norm(X(:,1));
Q(:,1)=X(:,1)/R(1,1); %��λ��
for k=2:m  
    for (i=1:k-1), R(i,k)=Q(:,i)'*X(:,k); end
    qt=X(:,k);
    for (i=1:k-1),qt=qt-R(i,k)*Q(:,i); end
    R(k,k)=norm(qt);  Q(:,k)=qt/R(k,k); %��λ��
end

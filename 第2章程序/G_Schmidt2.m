%改进Gram-Schmidt 正交化程序-G_Schmidt2.m
function [Q,R]=G_Schmidt2(X)
[n,m]=size(X); %矩阵的列
R=zeros(m);  Q=zeros(n,m);
R(1,1)=norm(X(:,1));
Q(:,1)=X(:,1)/R(1,1); %单位化
for k=2:m
    qt=X(:,k); 
    for i=1:k-1
        R(i,k)=qt'*Q(:,i);
        qt=qt-R(i,k)*Q(:,i); %对剩余向量进行修正
    end
    for i=1:k-1  %重正交化
        rt=qt'*Q(:,i); qt=qt-rt*Q(:,i); 
        R(i,k)=R(i,k)+rt;
    end
    R(k,k)=norm(qt);
    Q(:,k)=qt/R(k,k);%对本次得到的向量单位化
end


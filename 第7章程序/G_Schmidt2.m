%�Ľ�Gram-Schmidt ����������-G_Schmidt2.m
function [Q,R]=G_Schmidt2(X)
[n,m]=size(X); %�������
R=zeros(m);  Q=zeros(n,m);
R(1,1)=norm(X(:,1));
Q(:,1)=X(:,1)/R(1,1); %��λ��
for k=2:m
    qt=X(:,k); 
    for i=1:k-1
        R(i,k)=qt'*Q(:,i);
        qt=qt-R(i,k)*Q(:,i); %��ʣ��������������
    end
    for i=1:k-1  %��������
        rt=qt'*Q(:,i); qt=qt-rt*Q(:,i); 
        R(i,k)=R(i,k)+rt;
    end
    R(k,k)=norm(qt);
    Q(:,k)=qt/R(k,k);%�Ա��εõ���������λ��
end


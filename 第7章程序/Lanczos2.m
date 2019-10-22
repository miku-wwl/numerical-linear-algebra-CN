function [Q,T,betak] = Lanczos2(A,q1,k)
%给定对称矩阵A和单位长度向量q1和正整数k,
%计算一个长度为k的Lanczos分解:(重正交化)
%AQ_k=Q_kT_k+beta_kq_{k+1}e_k^T
q1=q1/norm(q1); Q(:,1)=q1;
w=A*Q(:,1); alpha(1)=Q(:,1)'*w;
w=w-alpha(1)*Q(:,1);
beta(1)=norm(w); Q(:,2)=w/beta(1);
for i=2:k
    w=A*Q(:,i); alpha(i)=Q(:,i)'*w;
    w=w-alpha(i)*Q(:,i)-beta(i-1)*Q(:,i-1);
    for j=1:i-1   %重正交化
        w=w-(Q(:,j)'*w)*Q(:,j);
    end
    beta(i)=norm(w);
    if (beta(i)<eps), 
        return;
    else
        Q(:,i+1)=w/beta(i);
    end
end
Q=Q(:,1:k); betak=beta(k);
T=diag(alpha)+diag(beta(1:k-1),-1)+diag(beta(1:k-1),1);



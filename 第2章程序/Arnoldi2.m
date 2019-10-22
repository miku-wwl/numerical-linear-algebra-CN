function [Q,H] = Arnoldi2(A,q1,k)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
q1=q1/norm(q1);  Q(:,1)=q1;
for j=1:k
    w=A*Q(:,j);
    for i=1:j
        H(i,j)=Q(:,i)'*w; w=w-H(i,j)*Q(:,i);
    end
    for i=1:j
        s=Q(:,i)'*w; H(i,j)=H(i,j)+s; w=w-s*Q(:,i);
    end
    H(j+1,j)=norm(w);
    if abs(H(j+1,j))<1.e-12
        return;
    else
        Q(:,j+1)=w/H(j+1,j);
    end
end



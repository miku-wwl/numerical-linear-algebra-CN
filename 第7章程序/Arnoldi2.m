function [V,H] = Arnoldi2(A,v,k)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
v=v/norm(v);  V(:,1)=v;
for j=1:k
    w=A*V(:,j);
    for i=1:j
        H(i,j)=V(:,i)'*w; w=w-H(i,j)*V(:,i);
    end
    for i=1:j
        s=V(:,i)'*w; H(i,j)=H(i,j)+s; w=w-s*V(:,i);
    end
    H(j+1,j)=norm(w);
    if abs(H(j+1,j))<1.e-12
        return;
    else
        V(:,j+1)=w/H(j+1,j);
    end
end



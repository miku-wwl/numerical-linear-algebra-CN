function [V,alpha,beta]=expand_lanczos(A,V,alpha,beta,m,k)
%给定对称矩阵A的长度为m的Lanczos分解,本算法将其扩展为一个长度为k的Lanczos分解
for i=m+1:k
    u=A*V(:,i)-beta(i-1)*V(:,i-1);
    alpha(i)=u'*V(:,i);
    u=u-alpha(i)*V(:,i);
    u=u-(u'*V(:,i-1))*V(:,i-1); %局部重正交化
    u=u-(u'*V(:,i))*V(:,i);
    u=u-V(:,1:i)*(V(:,1:i)'*u); %完全重正交化
    beta(i)=norm(u);
    if (beta(i)==0)
        return;
    else
        V(:,i+1)=u/beta(i);
        V=[V, V(:,i+1)];
    end
end


    
    



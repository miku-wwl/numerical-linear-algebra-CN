%不完全Cholesky分解程序
function [L] = michol(A)
%输入对称正定矩阵A,输出下三角矩阵L, 满足A=LL'+R
% n=8; e=ones(n,1); e1=ones(n-1,1); e2=ones(n-2,1);
%A=diag(8*e)-diag(2*e1,1)-diag(2*e1,-1)-diag(e2,2)-diag(e2,-2)；
n=size(A,1); L=zeros(n,n);
L(1,1)=sqrt(A(1,1)); L(2,1)=A(2,1)/L(1,1);
for (k=2:n)
    s=0;
    for (p=1:k-1)
        s=s+L(k,p)^2;
    end
    L(k,k)=sqrt(A(k,k)-s);
    for (i=k+1:n)
        if (A(i,k)~=0)
            s=0;
            for (p=1:k-1), s=s+L(i,p)*L(k,p); end
            L(i,k)=(A(i,k)-s)/L(k,k);
        end
    end
end
        
    

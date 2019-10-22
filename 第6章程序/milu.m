function [L,U,R]=milu(A,NF,omega)
%矩阵A的松弛不完全LU分解
%输入: 矩阵A, 非零模式NF, 松弛参数omega
%输出: 单位下三角阵L, 上三角阵U, 残差矩R, 满足A=LU+R
n=size(A,1);  R=zeros(n);
for k=1:n-1
    for i=k+1:n
        if (ismember(i+k*sqrt(-1),NF)==1)  %(i,k)属于集合NF
            A(i,k)=A(i,k)/A(k,k);
            for j=k+1:n
                A(i,j)=A(i,j)-A(i,k)*A(k,j);
            end
            for j=k+1:n
                if (ismember(i+j*sqrt(-1),NF)==0)  %(i,j)不属于集合NF
                    R(i,j)=R(i,j)+A(i,j);
                    R(i,i)=R(i,i)-omega*A(i,j);
                    A(i,i)=A(i,i)+omega*A(i,j);
                    A(i,j)=0;
                end
            end
        end
    end
end
L=tril(A,-1)+eye(n); U=triu(A);
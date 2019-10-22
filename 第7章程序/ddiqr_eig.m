function [iter,D]=ddiqr_eig(A,tol)
%本程序用双重位移隐式QR方法求实方阵的全部特征值.
%输入n阶实方阵A, 控制精度tol(默认是1e-5)
%输出迭代次数iter, A的全部特征值D
if nargin<2,tol=1e-5;end
n=size(A,1);
D=zeros(n,1); i=n; m=n; iter=0;  %初始化
[A]=hessenb(A); %化矩阵A为Hessenberg矩阵
while (m>0) 
    %用双重位移隐式QR方法进行迭代
    if m<=2
        la=eig(A(1:m,1:m)); 
        D(1:m)=la';  break;
    end
    iter=iter+1;
    A=ddiqr_tran(A,m);  %对上Hessenberg 矩阵做QR分解,并做正交相似变换
    for k=m-1:-1:1  %下面的程序段是判断是否终止
        if abs(A(k+1,k))<tol
            if m-k<=2
                la=eig(A(k+1:m,k+1:m));
                j=i-m+k+1; D(j:i)=la';
                i=j-1; m=k; break;
            end
        end
    end
end

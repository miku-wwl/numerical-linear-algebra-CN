%Givens-Householder方法程序-givens_householder.m
function [Lam,Q]=givens_householder(A,tol)
%本程序用Givens-Householder求实对称矩阵A的第m个特征值.
%输入参数:A为n阶实对称方阵, m特征值序号,tol为容许误差.
%输出参数:lambda第m个特征值,k为满足精度迭代次数
if nargin<4, max_it=500; end
if nargin<3, tol=1e-6; end
n=size(A,1);  k=0; Lam=zeros(n,1);
[T,Q]=mhessen(A); b=norm(T,inf); a=-b;
alpha=diag(T); beta=diag(T,1);
for i=1:n
    while(1)
        c=(a+b)/2;
        q=qkfun(c,alpha,beta);
        tt=find(q>=0); 
        sn=length(tt); 
        if(sn>=i), a=c; else, b=c; end
        if abs(b-a)<=tol,
            Lam(i)=c; break; 
        end
    end
end
   
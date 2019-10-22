%预处理共轭梯度法程序-pcg.m
function [x,iter,time,res,resvec]=pcg(A,b,x,M,max_it,tol)
%输入:系数矩阵A,右端向量b,初始向量x,预处
%理子M,容许误差tol最大迭代次数max_it
%输出:解向量x,迭代次数iter,CPU时间time,
%终止时相对残差模res相对残差模向量resvec
tic; r=b-A*x; z=M\r; p=z; 
rho=z'*r; mr=norm(r); iter=0;
while (iter<max_it)
    iter=iter+1;
    u=A*p; alpha=rho/(p'*u);
    x=x+alpha*p; r=r-alpha*u;
    z=M\r; rho1=z'*r;  
    beta=rho1/rho; p=z+beta*p;
    res=norm(r)/mr;
    resvec(iter)=res;
    if (res<tol), break; end
    rho=rho1;
end
time=toc;

%共轭梯度法程序-mcg.m
function [x,iter,time,res,resvec]=mcg(A,b,x,max_it,tol)
%输入:系数阵A,右端量b,初始值x,容许误差tol,最大迭代数max_it
%输出:解向量x,迭代次数iter,CPU时间time,
%终止时相对残差模res相对残差模向量resvec
tic; r=b-A*x; p=r; 
rho=r'*r; mr=sqrt(rho); iter=0;
while (iter<max_it)
    iter=iter+1;
    z=A*p; alpha=rho/(z'*p);
    x=x+alpha*p; r=r-alpha*z;
    rho1=r'*r; beta=rho1/rho;
    p=r+beta*p;
    res=sqrt(rho1)/mr;
    resvec(iter)=res;
    if (res<tol), break; end
    rho=rho1;
end
time=toc;

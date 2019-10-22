%BCG方法程序-bcg.m
function [x,k,time,res,resvec]=bcg(A,b,x,max_it,tol)
tic; n=length(b); q=zeros(n,1); qt=q; 
rho=1; r=b-A*x; mr=norm(r); k=0;
rt=ones(n,1);  %rt的选取会影响收敛速度 
%rt=r/mr; %rt=zeros(n,1); rt(1)=1;   
while (k<max_it)
    k=k+1;
    res=norm(r)/mr; resvec(k)=res;
    if (res<tol), break; end
    rho1=rt'*r; beta=rho1/rho;
    q=r+beta*q; qt=rt+beta*qt;
    Aq=A*q; sigma=qt'*Aq; 
    alpha=rho1/sigma;
    x=x+alpha*q; r=r-alpha*Aq;
    rt=rt-alpha*(A'*qt);
    rho=rho1;  
end
time=toc;

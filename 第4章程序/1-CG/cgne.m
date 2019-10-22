%CGNE·½·¨³ÌÐò-cgne.m
function [x,iter,time,res,resvec]=cgne(A,b,x,max_it,tol)
tic; r=b-A*x; p=A'*r;   
rho=r'*r; mr=sqrt(rho); iter=0;
while (iter<max_it)
    iter=iter+1; z=A*p;
    alpha=rho/(p'*p);
    x=x+alpha*p; r=r-alpha*z;
    rho1=r'*r; beta=rho1/rho;
    w=A'*r; p=w+beta*p; 
    res=norm(r)/mr;  
    resvec(iter)=res;
    if (res<tol), break; end
    rho=rho1;  
end
time=toc;

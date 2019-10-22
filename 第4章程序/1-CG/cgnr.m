%CGNR·½·¨³ÌÐò-cgnr.m
function [x,iter,time,res,resvec]=cgnr(A,b,x,max_it,tol)
tic; r=b-A*x; p=A'*r; u=A*p;
w=A'*r; rho=w'*w; mr=norm(r); iter=0;
while (iter<max_it)
    iter=iter+1;
    z=A'*p;  alpha=rho/(z'*z);
    x=x+alpha*p;  r=r-alpha*u;
    w=A'*r; rho1=w'*w; 
    beta=rho1/rho;
    p=w+beta*p; u=A*p;
    res=norm(r)/mr; 
    resvec(iter)=res;
    if (res<tol), break; end
    rho=rho1;
end
time=toc;

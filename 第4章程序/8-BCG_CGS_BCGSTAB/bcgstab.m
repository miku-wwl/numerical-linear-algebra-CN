%BCGSTAB·½·¨³ÌÐò-bcgstab.m
function [x,k,time,res,resvec]=bcgstab(A,b,x,max_it,tol)
tic; r=b-A*x;  p=r; mr=norm(r);
rt=ones(length(b),1); %rt=r; 
%rt=zeros(length(b),1);rt(1)=1;  
rho=rt'*r;  k=0;
while (k<=max_it)
    k=k+1; u=A*p; sigma=rt'*u; 
    alpha=rho/sigma; s=r-alpha*u;
    v=A*s; omega=(s'*v)/(v'*v);
    x=x+alpha*p+omega*s;
    r=s-omega*v; rho1=rt'*r; 
    beta=(alpha*rho1)/(omega*rho);
    p=r+beta*(p-omega*u);
    res=norm(r)/mr; resvec(k)=res;
    if (res<tol),  break; end
    rho=rho1;  
end
time=toc;

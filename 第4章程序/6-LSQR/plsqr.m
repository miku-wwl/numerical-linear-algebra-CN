%LSQR·½·¨³ÌÐò-mlsqr.m
function [x, k, time, res, resvec] = plsqr(A,b,x,M,max_it,tol)
tic; b=M\b; A=M\A; r=b-A*x; 
mr=norm(r); u=r/mr; v=A'*u; 
alpha=norm(v); v=v/alpha; z=v; 
zetat=mr; rhot=alpha; resvec(1)=1; k=0;
while (k<=max_it)
    k=k+1;
    u=A*v-alpha*u; beta=norm(u); u=u/beta;
    v=A'*u-beta*v; alpha=norm(v); v=v/alpha;
    rho=sqrt(rhot^2+beta^2);  c=rhot/rho; s=beta/rho;
    theta=s*alpha; rhot=c*alpha;
    zeta=c*zetat; zetat=-s*zetat;
    x=x+(zeta/rho)*z; r=b-A*x;    
    z=v-(theta/rho)*z;
    res=norm(r)/mr; resvec(k+1)=res;
    if (res<tol), break; end
end
time=toc;

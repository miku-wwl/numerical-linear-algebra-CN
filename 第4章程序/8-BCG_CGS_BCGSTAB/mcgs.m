function [x,k,time,res,resvec] = mcgs(A,b,x,max_it,tol)
%共轭梯度平方法(CGS)求解方程组Ax=b
tic; n=length(b); r=b-A*x;  
mr=norm(r);  resvec(1)=1;
%rt=ones(n,1);  %not converge
rt=zeros(n,1); rt(1)=1;  %rt=r; 
q=zeros(n,1); p=q; rho=1; k=0; 
while (k<max_it)
    rho1=rt'*r; beta=rho1/rho;
    u=r+beta*p; q=u+beta*(p+beta*q);
    v=A*q; sigma=rt'*v; alpha=rho1/sigma;
    p=u-alpha*v; z=alpha*(u+p);
    x=x+z;  r=r-A*z;
    res=norm(r)/mr; resvec(k+1)=res;
    if (res<tol)
        break;
    end
    rho=rho1; k=k+1;
end
time=toc;

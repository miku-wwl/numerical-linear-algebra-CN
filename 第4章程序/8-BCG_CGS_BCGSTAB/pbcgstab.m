function [x,k,time, res, resvec]=pbcgstab(A,b,x,M1,M2,max_it,tol)
%共轭梯度法求解对称正定方程组Ax=b
tic; r=b-A*x;  p=r; mr=norm(r);
rt=r; %rt=ones(length(b),1);
%rt=zeros(length(b),1); rt(1)=1; 
rho=rt'*r;  k=1; 
while (k<=max_it)
    y=M2\(M1\p); u=A*y;  sigma=rt'*u; 
    alpha=rho/sigma; s=r-alpha*u;
    z=M2\(M1\s); v=A*z;  
    xi=M1\v; eta=M1\s;
    omega=(xi'*eta)/(xi'*xi);
    x=x+alpha*y+omega*z;
    r=s-omega*v;
    rho1=rt'*r; beta=(alpha*rho1)/(omega*rho);
    p=r+beta*(p-omega*u);
    res=norm(r)/mr;
    resvec(k)=res;
    if (res<tol),  break;  end
    k=k+1;
    rho=rho1;
end
time=toc;

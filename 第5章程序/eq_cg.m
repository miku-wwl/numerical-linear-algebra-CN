function [x, k, time, res, resvec] = eq_cg(A,b,x,max_it,tol)
%共轭梯度法求解对称正定方程组Ax=b
if nargin<5, max_it=1000; end
if nargin<4, tol=1.e-6; end
if nargin<3, x=zeros(size(b)); end
r0=b-A*x;  p=r0; 
rho=r0'*r0;  mr0=sqrt(rho);  k=0;
tic;
while (k<max_it)
    k=k+1;
    z=A*p; nu=rho/(z'*p);
    x=x+nu*p;  r=r0-nu*z;
    rho1=r'*r; mu=rho1/rho;
    p1=r+mu*p;
    res=sqrt(rho1)/mr0;
    resvec(k)=res;
    if (res<tol)
        break;
    end
    r0=r; p=p1; rho=rho1;
end
time=toc;

%递推形式的极小残量法-minresr.m
function [x,iter,time,res,resvec]=minresr(A,b,x,max_it,tol)
%极小残量法求解对称不定方程组Ax=b
tic; r=b-A*x;  p=r; mr=norm(r); 
u=A*p; z=A*r; iter=1;
while (iter<max_it)
    alpha=(r'*u)/(u'*u);
    x=x+alpha*p;  r1=r-alpha*u; 
    z1=A*r1;  beta=(z1'*r1)/(z'*r);
    p=r1+beta*p;   u=z1+beta*u;
    res=norm(r1)/mr; resvec(iter)=res;
    if (res<tol),   break;  end
    r=r1; z=z1; iter=iter+1; 
end
time=toc;


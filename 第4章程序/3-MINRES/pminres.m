function [x,k,time,res,resvec] = pminres(A,b,x,M,max_it,tol)
%预处理极小残量法求解对称不定方程组Ax=b
n=length(b);
r=b-A*x;  z0=M\r; p=z0; mr=norm(r);
u0=A*p; v0=u0; w0=M\u0; k=1;
tic;
while (k<max_it)
    alpha=(z0'*u0)/(w0'*u0); x=x+alpha*p;
    r=r-alpha*u0;   z=M\r;   v=A*z;
    beta=(v'*z)/(v0'*z0);    p=z+beta*p;
    u=v+beta*u0;    w=M\u;
    res=norm(r)/mr;   resvec(k)=res;
    if (res<tol),   break;  end
    z0=z; v0=v;
    u0=u; w0=w;  
    k=k+1;
end
time=toc;


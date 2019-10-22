%PMINRES方法程序-pminres.m
function [x,iter,time,res,resvec] = pminres1(A,b,x,M1,M2,max_it,tol)
%PMINRES方法求解对称不定方程组Ax=b, 预处理子M=M1*M2
tic; n=length(b);
r=b-A*x; z=M2\(M1\r); p=z; mr=norm(r);
u=A*p; v=u; w=M2\(M1\u); iter=1;
while (iter<max_it)
    alpha=(z'*u)/(w'*u); x=x+alpha*p;
    r=r-alpha*u; z1=M2\(M1\r); v1=A*z1;
    beta=(v1'*z1)/(v'*z); p=z1+beta*p;
    u=v1+beta*u; w=M2\(M1\u);
    res=norm(r)/mr; resvec(iter)=res;
    if (res<tol),  break;  end
    z=z1; v=v1; iter=iter+1;
end
time=toc;


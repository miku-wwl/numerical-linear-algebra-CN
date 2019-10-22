%MINRES方法的程序-mminres.m
function [x,iter,time,res,resvec]=mminres(A,b,x,max_it,tol)
%极小残量法求解对称不定方程组Ax=b
tic; r=b-A*x; b=norm(r); v=r/b; bt=b;
z=A*v; a=v'*z; u=z-a*v; b1=norm(u); 
if (b1~=0),  v1=u/b1; end
c0=1; s0=0; p0=zeros(length(b),1);
[c,s,gama]=givens(a,b1); %Givens变换
p=v/gama; rho=-b*s; tau=b*c;
x=x+tau*p;  iter=1;  
while (iter<max_it)
    res=abs(rho)/bt; resvec(iter)=res;
    if (res<tol),   break;  end
    z=A*v1; a=v1'*z; u=z-a*v1-b1*v; b2=norm(u);
    if (b2~=0),  v2=u/b2;   end
    epsi=s0*b1; bh1=c0*b1;
    dta=c*bh1+s*a; ah=-s*bh1+c*a;
    [c1,s1,gama]=givens(ah,b2); %Givens变换
    tau=rho*c1; rho=-rho*s1;
    p1=(v1-epsi*p0-dta*p)/gama;
    x=x+tau*p1; b1=b2; iter=iter+1; 
    v=v1; v1=v2; p0=p; p=p1;  
    c0=c; c=c1; s0=s; s=s1; 
end
time=toc;


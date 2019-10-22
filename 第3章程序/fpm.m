%clc
function [x,iter,time,res]=fpm(A,b,x0, tol)
n=length(b); In=eye(n); 
B=In-A;  iter=0; x=x0;
res=norm(b-A*x)/norm(b);
tic
while (iter<1000)
    iter=iter+1;
    x=B*x+b;
    res=norm(b-A*x)/norm(b);   
    if res<tol
       break;
    end
end
time=toc;


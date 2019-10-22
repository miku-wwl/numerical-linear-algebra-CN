function [x,fval,k,time]=nls_cg(A,b,x,tol,max_it)
%CG算法求解法方程A'Ax=A'b
tic; n=size(A,2); r=b-A*x; s=A'*r; 
p=s; gama=s'*s; nr=norm(s); k=0;
while (k<=max_it)
    k=k+1;
    q=A*p; alpha=gama/(q'*q);
    x=x+alpha*p; 
    r=r-alpha*q; s=A'*r;
    if (norm(s)/nr<tol), break; end
    gama1=s'*s; beta=gama1/gama;
    p=s+beta*p;
    gama=gama1; 
end
fval=norm(r);
time=toc;


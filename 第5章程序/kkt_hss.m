function [x,fval,k,time]=kkt_hss(A,b,alpha,x,tol,max_it)
%HSS算法求解KKT方程r=b-Ax, A'r=0
if nargin<6, max_it=1000; end
if nargin<5, tol=1.e-6; end
if nargin<4, x=zeros(size(A,2),1); end
tic; m=size(A,1); r=b-A*x; 
s=A'*r; nr=norm(s); k=0;
B=alpha*eye(m)+A*A'/alpha;
while (k<=max_it)
    k=k+1;
    r=(alpha*r-A*x+b)/(alpha+1);
    x=x+s/alpha; 
    c=(alpha-1)*r-A*x+b; 
    r=B\c;  %精确求解子问题
    s=A'*r;  x=x+s/alpha;
    if (norm(s)/nr<tol), break; end
end
fval=norm(r);
time=toc;


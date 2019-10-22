function [x,fval,k,time]=kkt_sor(A,b,x,ome, tol,max_it)
%Jacobi算法求解KKT方程
tic; [m,n]=size(A); r=b-A*x; v=r(1:n); w=r(n+1:m);
A1=A(1:n,:); A2=A(n+1:m,:); A3=A2/A1;
b1=b(1:n); b2=b(n+1:m);  nr=norm(A'*r); k=0; 
tau=ome*(1-ome);
while (k<=max_it)
    k=k+1;
    x=(1-ome)*x+ome*A1\(b1-v); 
    w=(1-ome)*w-tau*A2*x-ome^2*A3*(b1-v)+ome*b2; 
    v=(1-ome)*v+tau*(1-ome)*A3'*A2*x-tau*A3'*w-ome^2*A3'*(b2-ome*A3*(b1-v));
    r=b-A*x;  s=A'*r;
    if (norm(s)/nr<tol)
        break;
    end 
end
fval=norm(r); time=toc;


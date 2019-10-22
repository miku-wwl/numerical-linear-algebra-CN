function [x,fval,k,time]=kkt_jacobi(A,b,x,tol,max_it)
%Jacobi算法求解KKT方程
tic; [m,n]=size(A); r=b-A*x; v=r(1:n); w=r(n+1:m);
A1=A(1:n,:); A2=A(n+1:m,:); %A3=A2/A1;
b1=b(1:n); b2=b(n+1:m);  nr=norm(A'*r); k=0;
B=[A1,zeros(n,m-n),eye(n); A2, eye(m-n), zeros(m-n,n); zeros(n), A2',A1']; 
while (k<=max_it)
    k=k+1;
    x=A1\(b1-v); w=b2-A2*x; 
    v=-A1'\(A2'*w);
    r=b-A*x;  s=A'*r;
    if (norm(s)/nr<tol)
        break;
    end 
end
fval=norm(r); time=toc;


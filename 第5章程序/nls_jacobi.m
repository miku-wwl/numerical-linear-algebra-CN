function [x,fval,k,time]=nls_jacobi(A,b,x,tol,max_it)
%Jacobi算法求解法方程A'Ax=A'b
tic; n=size(A,2); r=b-A*x; 
d=zeros(n,1); nr=norm(A'*r); k=0;
for i=1:n
    d(i)=A(:,i)'*A(:,i);
end
%D=diag(diag(d)) ;
while (k<=max_it)
    k=k+1;
    %x=D\((D-B)*x+c);
    for i=1:n
       x(i)=x(i)+A(:,i)'*r/d(i);
    end
    %x=x(:),  
    r=b-A*x; s=A'*r;
    if (norm(s)/nr<tol)
        break;
    end 
end
fval=norm(r); time=toc;


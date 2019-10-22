function [x,fval,k,time]=nls_sor(A,b,w,x,tol,max_it)
%SOR迭代法求解法方程A'Ax=A'b
tic; n=size(A,2); r=b-A*x; d=zeros(n,1);
nr=norm(A'*r); k=0; E=eye(n);
for i=1:n
    d(i)=A(:,i)'*A(:,i);
end
D=diag(d) ; Z=[x]; R=[r];
while (k<=max_it)
    k=k+1;
    for i=1:n
        delta=w*A(:,i)'*R(:,i)/d(i);
        z=Z(:,i)+delta*E(:,i);
        r=R(:,i)-delta*A(:,i);
        %x(i)=x(i)+A(:,i)'*r/d(i);
        Z=[Z,z]; R=[R,r];
    end
    x=Z(:,n+1);
    r=b-A*x;  Z=[x]; R=[r];
    if (norm(A'*r)/nr<tol)
        break;
    end
end
fval=norm(r);
time=toc;


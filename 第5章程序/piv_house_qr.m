function [x, fval, P]=piv_house_qr(A,b)
%秩亏最小二乘问题的列主元QR分解法
[m,n]=size(A); %Q=eye(m); 
A1=A; b1=b;
P=eye(n);
for j=1:n
    c(j)=A(1:m,j)'*A(1:m,j);  
end
[cr,r]=max(c);
for k=1:n
    if cr<=0, break; end  
    c([k r])=c([r k]); 
    P(:,[k r])=P(:,[r k]);   
    A(1:m, [k r])= A(1:m, [r k]);
    [v,beta]=house(A(k:m,k));
    H=eye(m-k+1)-beta*v*v';
    A(k:m,k:n)=H*A(k:m,k:n);
    b(k:m)=H*b(k:m);
    %Q=Q*blkdiag(eye(k-1),(eye(m-k+1)-beta*v*v'));
    %A(k+1:m,k)=v(2:m-k+1);
    for j=k+1:n
        c(j)=c(j)-A(k,j)^2;
    end
    [cr,r]=max(c(k+1:n));
    r=r+k;
end
for i=1:m
    if sum(abs(A(i,:)))<1.0e-10
        rank=i-1; break;
    end
end
R11=A(1:rank,1:rank); c=b(1:rank);
x=P*[R11\c; zeros(rank,1)];
% S=A(1:rank,:)*P';  x=S\c;
fval=norm(A1*x-b1);


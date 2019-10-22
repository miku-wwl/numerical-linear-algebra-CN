%Àý5.3-ex53.m
A=[1 2 3 4; 1 4 5 6; 1 5 6 7; 1 8 9 10; 1 11 12 13];
b=[11 13 15 18 20]';
[m,n]=size(A);
[U,S,V]=svd(A);
% x=(V*pinv(S)*U')*b
for i=1:n
    if abs(S(i,i))<1.0e-6
        r=i-1; break;
    end
end    %r=rank(S);
x=zeros(n,1);
for i=1:r
    x=x+(U(:,i)'*b/S(i,i))*V(:,i);
end
x
fval=norm(A*x-b)
% xb=A\b
% fval=norm(A*xb-b)
    
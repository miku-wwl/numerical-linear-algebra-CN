%Àý7.10 
clear all
n=100; A=zeros(n,n);
k=50; I=eye(k); q=rand(n,1);
for i=1:n
    A(i,i)=i;
end
for i=1:n-1,
    for j=i+1:n
        A(i,j)=i;  
    end
end
for i=2:n
    for j=1:i-1
        A(i,j)=-j;
    end
end
tic, [Q1,H1] = Arnoldi(A,q,k);  toc
tic, [Q2,H2] = Arnoldi2(A,q,k);  toc
err1=norm(Q1(:,1:k)'*Q1(:,1:k)-I,'fro')
err2=norm(Q2(:,1:k)'*Q2(:,1:k)-I,'fro')
d=eig(A);
[s,t]=max(abs(d));
d=d(t)



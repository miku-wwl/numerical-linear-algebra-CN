%Àý 
clear all
n=1000; A=zeros(n,n);
k=50; I=eye(k); q=rand(n,1);
for i=1:n,
    for j=1:i
        A(i,j)=j;  A(j,i)=A(i,j);
    end
end
tic,  [Q1,T1,betak1] = Lanczos(A,q,k);  toc
tic,  [Q2,T2,betak2] = Lanczos2(A,q,k);  toc
err1=norm(Q1'*Q1-I,'fro')
err2=norm(Q2'*Q2-I,'fro')



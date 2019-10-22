%Àý2.7-ex27
clear all
n=1000; A=zeros(n,n);
k=50; I=eye(k); v=rand(n,1);
for i=1:n,
    for j=1:i
        A(i,j)=j;  A(j,i)=A(i,j);
    end
end
tic,  [V1,T1,betak1] = Lanczos(A,v,k);  toc
tic,  [V2,T2,betak2] = Lanczos2(A,v,k);  toc
err1=norm(V1'*V1-I,'fro')
err2=norm(V2'*V2-I,'fro')



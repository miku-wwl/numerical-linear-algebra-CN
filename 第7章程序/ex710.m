%Àı7.10 
n=12; tol=1.e-6; max_it=500;
A=zeros(n); lambda=zeros(n,1); iter=zeros(n,1);
for i=1:n,
    for j=1:i
        A(i,j)=j;  A(j,i)=A(i,j);
    end
end
[T,Q]=mhessen(A);
for m=1:n
    [lambda(m),iter(m)]=givens_househ(T,m,tol,max_it);
end
D=[sort(lambda),    eig(A)];
disp('    ÌØÕ÷Öµ,   eigÃüÁî')
disp(D)
%iter
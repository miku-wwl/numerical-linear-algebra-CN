%A和M\A的谱分布
n=100; tol=1e-21;
e=[1:n]'; e1=ones(n-1,1);
A=diag(e)-diag(e1,1)-diag(e1,-1);
M=diag(diag(A));
A1=M\A;
ev1=eig(A); ev2=eig(A1);
t=1:2:n;
v1=ev1(t); v2=ev2(t);
plot(t, v1,'r*'); hold on
plot(t,v2,'bx');
legend('A 的谱分布', 'M^{-1}A 的谱分布');
axis([1,100,-1, 120])


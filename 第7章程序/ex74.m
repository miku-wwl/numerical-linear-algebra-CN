%例7.4
clear all
n=100; e=[1:n]'; e1=ones(n-1,1);
A=diag(e)-diag(e1,1)-diag(e1,-1);
[d,V]=Jacobi_eig(A);
lam=eig(A);
err=norm(d-lam,inf)
tt=1:2:n;
plot(tt,lam(tt),'*-k');
axis([0,101,0,101]);
ylabel('特征值')
title('矩阵的特征值分布')



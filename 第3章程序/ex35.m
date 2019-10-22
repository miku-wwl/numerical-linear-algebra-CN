%例3.5-ex3.5
n=2^(14)-1; e=ones(n,1); 
A=spdiags([-e, 4*e, -e],[-1 0 1],n,n);
b=2*e; b(1)=3; b(end)=3;
tol=1.e-10; x=zeros(n,1);
lamd=0.5*cos(pi/(n+1)); m=10;
[x1,iter1,res1,t1]=eig_extr(A,b,x,m,lamd,tol);
[x2,iter2,res2,t2]=mjacobi(A,b,x,tol);
%格式化显示
fid = 1;
fprintf(fid, '     算  法             Iter      CPU            RES\n' );     fprintf('\n');
fprintf(fid, '特征值外推法,      %4i    %4.4f     %11.4e\n', iter1, t1,res1); fprintf('\n');
fprintf(fid, 'Jacobi迭代法,     %4i    %4.4f     %11.4e\n', iter2, t2, res2); fprintf('\n');


%例3.6-ex36.m
n=2^(14)-1; e=ones(n,1); 
A=spdiags([-e, 4*e, -e],[-1 0 1],n,n);
b=2*e; b(1)=3; b(end)=3;
tol=1.e-10; x0=zeros(n,1);
h=1.0/(n+1); b1=0.5*cos(pi*h); a1=0.5*cos(n*pi*h);
%a1=-0.9; b1=0.9;  %a1=-1; b1=1;  
[x1,iter1,res1,t1]=mjacobi(A,b,x0,tol);
[x2,iter2,res2,t2]=mcheby(A,b,x0,a1,b1,tol);
%格式化显示
fid = 1;
fprintf(fid, '     算  法             Iter      CPU            RES\n' );     fprintf('\n');
fprintf(fid, '  Jacobi迭代法,     %4i    %4.4f     %11.4e\n', iter1,t1,res1); fprintf('\n');
fprintf(fid, 'Chebyshev加速,   %4i    %4.4f     %11.4e\n', iter2, t2,res2); fprintf('\n');


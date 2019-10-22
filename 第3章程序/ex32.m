%例3.2-ex32
n=4095; w=1.1;  %w=1.6;
tol=1.e-10; x0=zeros(n,1);
e=ones(n,1); 
A=spdiags([-e, 4*e, -e],[-1 0 1],n,n);
b=2*e; b(1)=3; b(end)=3;
[x1, iter1,err1, t1] =mjacobi(A, b, x0, tol);
[x2, iter2,err2, t2] = mseidel(A, b, x0, tol);
[x3, iter3,err3,t3] = msor(A, b, w, x0, tol);
[x4, iter4,err4, t4] = mssor( A, b, w, x0, tol);
%格式化显示
fid = 1;
fprintf(fid, '算法     Iter     CPU   ERR\n' );     fprintf('\n');
fprintf(fid, 'Jacobi, %4i   %8.4f     %11.4e\n', iter1,t1,err1); fprintf('\n');
fprintf(fid, 'G-S,    %4i    %8.4f     %11.4e\n', iter2, t2,err2); fprintf('\n');
fprintf(fid, 'SOR,   %4i    %8.4f     %11.4e\n', iter3,t3, err3); fprintf('\n');
fprintf(fid, 'SSOR, %4i    %8.4f     %11.4e\n', iter4, t4, err4); fprintf('\n');


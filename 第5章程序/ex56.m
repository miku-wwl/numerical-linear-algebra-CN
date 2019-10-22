%例5.6-ex56.m
clear all
A=[1,   23.73,   5.49,   1.21; 1,  22.34,   4.32,   1.35;  ...     
      1,  28.84,   5.04,   1.92;  1,  27.67,   4.72,   1.49; ...
      1,  20.83,   5.35,   1.56;  1,  22.27,   4.27,   1.50; ...
      1,  27.57,   5.25,   1.85;  1,  28.01,   4.62,   1.51];
b=[15.02, 12.62, 14.86, 13.98, 15.91, 12.47, 15.80, 14.32]';
% A=[1 2 3 4; 1 4 5 6; 1 5 6 7; 1 8 9 10; 1 11 12 13];
% b=[11 13 15 18 20]';
x=zeros(size(A,2),1); tol=1e-6; max_it=3000; 
w=1.06;
%[x1,fval1,k1,t1]=nls_jacobi(A,b,x,tol,max_it);
[x2,fval2,k2,t2]=nls_seidel(A,b,x,tol,max_it); 
[x3,fval3,k3,t3]=nls_sor(A,b,w,x,tol,max_it);
fprintf('\n'); fid = 1;
fprintf(fid, 'SOR迭代法     迭代次数    极小值     CPU时间   ' );     fprintf('\n');
fprintf(fid, 'omega=1       %4i      %8.4f     %8.4f\n', k2, fval2, t2);   
fprintf(fid, 'ome=1.06      %4i      %8.4f     %8.4f\n',  k3, fval3, t3);   fprintf('\n');
x3
%例5.4-ex54.m
clear all
A=[4,-1,0,0; 1, 4, -1, 0;  0,1, 4,-1;...
      0, 0,1,4; 1,3,2,1; 0,2,0,3; 8, 2, 3, 1];
b=[9,12,11,13,17,15,19]';
x=zeros(size(A,2),1); tol=1e-10; max_it=2000; 
[x1,fval1,k1,t1]=nls_jacobi(A,b,x,tol,max_it);
fprintf('\n'); fid = 1;
fprintf(fid, 'Jacobi    迭代次数    极小值     CPU时间   ' );     fprintf('\n');
fprintf(fid, '               %4i      %8.4f     %8.4f\n', k1, fval1, t1);  
x1
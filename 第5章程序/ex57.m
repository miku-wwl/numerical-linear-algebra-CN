%例5.7-ex57.m
clear all
A=[1,   23.73,   5.49,   1.21; 1,  22.34,   4.32,   1.35;  ...     
      1,  28.84,   5.04,   1.92;  1,  27.67,   4.72,   1.49; ...
      1,  20.83,   5.35,   1.56;  1,  22.27,   4.27,   1.50; ...
      1,  27.57,   5.25,   1.85;  1,  28.01,   4.62,   1.51];
b=[15.02, 12.62, 14.86, 13.98, 15.91, 12.47, 15.80, 14.32]';
x=zeros(size(A,2),1); tol=1e-10; max_it=5000;   
[x1,fval1,k1,time1]=nls_sor(A,b,1.06,x,tol,max_it); 
[x,fval,k,time]=nls_cg(A,b,x,tol,max_it);
fprintf('\n'); fid = 1;
fprintf(fid, ' 算法   迭代次数    极小值     CPU时间   ' );     fprintf('\n');
fprintf(fid, ' SOR    %4i      %8.4f     %8.4f\n', k1, fval1, time1);   
fprintf(fid, ' CG      %4i        %8.4f     %8.4f\n', k, fval, time);    fprintf('\n');
disp('    x1=        x=');
disp([x1,x]);


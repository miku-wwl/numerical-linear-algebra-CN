%例5.8-ex58.m
clear all
A=[1,   23.73,   5.49,   1.21; 1,  22.34,   4.32,   1.35;  ...     
      1,  28.84,   5.04,   1.92;  1,  27.67,   4.72,   1.49; ...
      1,  20.83,   5.35,   1.56;  1,  22.27,   4.27,   1.50; ...
      1,  27.57,   5.25,   1.85;  1,  28.01,   4.62,   1.51];
b=[15.02, 12.62, 14.86, 13.98, 15.91, 12.47, 15.80, 14.32]';
n=size(A,2);
x=zeros(n,1); tol=1e-10; max_it=10000;  alpha=6.0;
%[U,S,V]=svd(A);  alpha=sqrt(S(1,1)*S(n,n));
[x1,fval1,k1,time1]=kkt_ihss(A,b,alpha,x,tol,max_it);
[x2,fval2,k2,time2]=kkt_hss(A,b,alpha,x,tol,max_it);
[x3,fval3,k3,time3]=nls_cg(A,b,x,tol,max_it);
fprintf('\n'); fid = 1;
fprintf(fid, ' 算法   迭代次数    极小值     CPU时间   ' );     fprintf('\n');
fprintf(fid, 'IHSS    %4i      %8.4f     %8.4f\n', k1, fval1, time1); 
fprintf(fid, ' HSS    %4i      %8.4f     %8.4f\n', k2, fval2, time2); 
fprintf(fid, ' CG     %4i        %8.4f     %8.4f\n', k3, fval3, time3);    fprintf('\n');
disp('    x1=        x2=      x3 ');
disp([x1,x2,x3]);



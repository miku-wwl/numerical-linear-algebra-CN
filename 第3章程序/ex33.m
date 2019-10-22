%例3.3-ex3.3
m=2;  s=2; b=[1 1 1]'; tol=1.e-5; 
x0=[0 1 2]';
A=[1.2 0.3 0.4; 0.4 1.2 0.3; 0.3 0.4 1.2];
[x1,iter1,t1,res1]=fpm(A,b,x0, tol);
[x2,iter2,t2,res2]=wcm(A,b,x0,m,s,tol);
%格式化显示
fid = 1;
fprintf(fid, '    算法          Iter     CPU     RES\n' );     fprintf('\n');
fprintf(fid, 'Richardson,   %4i    %4.4f     %11.4e\n', iter1,t1,res1); fprintf('\n');
fprintf(fid, '整体校正法,   %4i    %4.4f     %11.4e\n', iter2, t2,res2); fprintf('\n');
x1


%例3.4-ex3.4
n=1023; e=ones(n,1); 
A=spdiags([-e, 4*e, -e],[-1 0 1],n,n);
b=2*e; b(1)=3; b(end)=3;
m=10;  s=10;  tol=1.e-10; 
x0=zeros(n,1);
[x1,iter1,t1,res1]=fpm(A,b,x0, tol);
[x2,iter2,t2,res2]=wcm(A,b,x0,m,s,tol);
[x3,iter3,t3,res3]=wcmj(A,b,x0,m,s,tol);
%格式化显示
fid = 1;
fprintf(fid, '    算法          Iter     CPU            RES\n' );     fprintf('\n');
fprintf(fid, 'Richardson,   %4i    %4.4f     %11.4e\n', iter1,t1,res1); fprintf('\n');
fprintf(fid, '整体校正法1,   %4i    %4.4f     %11.4e\n', iter2, t2,res2); fprintf('\n');
fprintf(fid, '整体校正法2,   %4i    %4.4f     %11.4e\n', iter3, t3,res3); fprintf('\n');




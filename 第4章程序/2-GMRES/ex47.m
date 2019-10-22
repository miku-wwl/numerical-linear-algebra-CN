%例4.7-ex47.m
clear all;  
n=1000; e=[1:n]'; e1=ones(n-1,1); 
A=diag(e)+diag(e1,-1)-diag(e1,1); A(1,n)=n; A(n,1)=-n;
%load SiNa; A=Problem.A;  n=size(A,1);  
M=diag(diag(A)); %预处理子
xs=ones(n,1); b=A*xs; x=zeros(n,1); tol=1e-10;
max_it=200;    rm=norm(b);
[x1,k1,j1,ti1,res1,resvec1,flag1]=pgmresm(A,b,x,M,3,max_it,tol);
tic,[x2, flag2, res2, k2, resvec2] = gmres(A, b, 3, tol, max_it, M); ti2=toc;  %系统自带的GMRES函数
[x3,k3,j3,ti3,res3,resvec3,flag3]=pgmresm(A,b,x,M,6,max_it,tol);
tic,[x4, flag4, res4, k4, resvec4] = gmres(A, b, 6, tol, max_it, M); ti4=toc;  %系统自带的GMRES函数
[x5,k5,j5,ti5,res5,resvec5,flag5]=pgmresm(A,b,x,M,9,max_it,tol);
tic,[x6, flag6, res6, k6, resvec6] = gmres(A, b, 9, tol, max_it, M); ti6=toc;  %系统自带的GMRES函数
[x7,k7,j7,ti7,res7,resvec7,flag7]=pgmresm(A,b,x,M,12,max_it,tol);
tic,[x8, flag8, res8, k8, resvec8] = gmres(A, b, 12, tol, max_it, M); ti8=toc;  %系统自带的GMRES函数
[x0, k0,  ti0, res0, resvec0, flag0]=pgmres (A, b, x, M, max_it, tol );       %PGMRES,
%显示误差
err0=norm(xs-x0); err1=norm(xs-x1); err2=norm(xs-x2); err3=norm(xs-x3); err4=norm(xs-x4);
err5=norm(xs-x5); err6=norm(xs-x6); err7=norm(xs-x7); err8=norm(xs-x8);
%格式化显示
fprintf('\n');  fid = 1;
fprintf(fid, ' m  Outer   Iner  Iter   flag      Time          Res          Err\n' );     fprintf('\n');
fprintf(fid, ' 3,  %4i   %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k1, j1, (k1-1)*3+j1, flag1, ti1, res1, err1); fprintf('\n');
fprintf(fid, ' 3,  %4i   %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k2(1),k2(2), (k2(1)-1)*3+k2(2), flag2, ti2, res2, err2); fprintf('\n');
fprintf(fid, ' 6,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k3, j3, (k3-1)*6+j3, flag3, ti3, res3, err3); fprintf('\n');
fprintf(fid, ' 6,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k4(1),k4(2), (k4(1)-1)*6+k2(2), flag4, ti4, res4, err4); fprintf('\n');
fprintf(fid, ' 9,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k5, j5, (k5-1)*9+j5, flag5, ti5, res5, err5); fprintf('\n');
fprintf(fid, ' 9,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k6(1),k6(2), (k6(1)-1)*9+k6(2), flag6, ti6, res6, err6); fprintf('\n');
fprintf(fid, '12, %4i   %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k7, j7, (k7-1)*12+j7, flag7, ti7, res7, err7); fprintf('\n');
fprintf(fid, '12, %4i   %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k8(1),k8(2), (k8(1)-1)*12+k8(2), flag8, ti8, res8, err8); fprintf('\n');
fprintf(fid, 'PG %4i   %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', 1,k0, k0, flag0, ti0, res0, err0); fprintf('\n');
%可视化
t1=1:length(resvec1);  t3=1:length(resvec3);
t5=1:length(resvec5);  t7=1:length(resvec7);
t0=1:length(resvec0);
semilogy(t1,resvec1(t1),'.-k');  hold on
semilogy(t3,resvec3(t3),'*-k');   
semilogy(t5,resvec5(t5),'x-k');  
semilogy(t7,resvec7(t7),'+-k');  
semilogy(t0,resvec0(t0),'o-k') ;hold off
legend('PGMRES(3)','PGMRES(6)','PGMRES(9)','PGMRES(12)','PGMRES')
xlabel('迭代次数'); ylabel('相对残差范数');
%axis([1,k0+10,0.5e-11, 1])


 

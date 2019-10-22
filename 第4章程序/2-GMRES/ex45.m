%例4.5-ex45.m
clear all;  
n=1000; e=[1:n]'; e1=ones(n-1,1); 
A=diag(e)+diag(e1,-1)-diag(e1,1); A(1,n)=n; A(n,1)=-n;
%load SiNa; A=Problem.A;  n=size(A,1);  
xs=ones(n,1); b=A*xs; x=zeros(n,1); tol=1e-10;
max_it=200;    rm=norm(b);
[x1,k1,j1,ti1,res1,resvec1,flag1]=gmresm(A,b,x,10,max_it,tol);
tic,[x2, flag2, res2, k2, resvec2] = gmres(A, b, 10, tol, max_it); ti2=toc;  %系统自带的GMRES函数
[x3,k3,j3,ti3,res3,resvec3,flag3]=gmresm(A,b,x,20,max_it,tol);
tic,[x4, flag4, res4, k4, resvec4] = gmres(A, b, 20, tol, max_it); ti4=toc;  %系统自带的GMRES函数
[x5,k5,j5,ti5,res5,resvec5,flag5]=gmresm(A,b,x,30,max_it,tol);
tic,[x6, flag6, res6, k6, resvec6] = gmres(A, b, 30, tol, max_it); ti6=toc;  %系统自带的GMRES函数
[x7,k7,j7,ti7,res7,resvec7,flag7]=gmresm(A,b,x,40,max_it,tol);
tic,[x8, flag8, res8, k8, resvec8] = gmres(A, b, 40, tol, max_it); ti8=toc;  %系统自带的GMRES函数
[x9,k9,j9,ti9,res9,resvec9,flag9]=gmresm(A,b,x,50,max_it,tol);
tic,[x10,flag10,res10,k10,resvec10] = gmres(A, b, 50, tol, max_it); ti10=toc;  %系统自带的GMRES函数
[x11,k11,j11,ti11,res11,resvec11,flag11]=gmresm(A,b,x,60,max_it,tol);
tic,[x12,flag12,res12,k12,resvec12] = gmres(A, b, 60, tol, max_it); ti12=toc;  %系统自带的GMRES函数
[x0, k0,  ti0, res0, resvec0, flag0] = mgmres (A, b, x, max_it, tol );       %GMRES,
%显示误差
err0=norm(xs-x0); err1=norm(xs-x1); err2=norm(xs-x2); err3=norm(xs-x3); err4=norm(xs-x4);
err5=norm(xs-x5); err6=norm(xs-x6); err7=norm(xs-x7); err8=norm(xs-x8);
err9=norm(xs-x9); err10=norm(xs-x10); err11=norm(xs-x11); err12=norm(xs-x12);
%格式化显示
fprintf('\n');  fid = 1;
fprintf(fid, ' m  Outer   Iner  Iter   flag      Time          Res          Err\n' );     fprintf('\n');
fprintf(fid, '10,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k1, j1, (k1-1)*10+j1, flag1, ti1, res1, err1); fprintf('\n');
fprintf(fid, '10,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k2(1),k2(2), (k2(1)-1)*10+k2(2), flag2, ti2, res2, err2); fprintf('\n');
fprintf(fid, '20,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k3, j3, (k3-1)*20+j3, flag3, ti3, res3, err3); fprintf('\n');
fprintf(fid, '20,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k4(1),k4(2), (k4(1)-1)*20+k2(2), flag4, ti4, res4, err4); fprintf('\n');
fprintf(fid, '30,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k5, j5, (k5-1)*30+j5, flag5, ti5, res5, err5); fprintf('\n');
fprintf(fid, '30,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k6(1),k6(2), (k6(1)-1)*30+k6(2), flag6, ti6, res6, err6); fprintf('\n');
fprintf(fid, '40,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k7, j7, (k7-1)*40+j7, flag7, ti7, res7, err7); fprintf('\n');
fprintf(fid, '40,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k8(1),k8(2), (k8(1)-1)*40+k8(2), flag8, ti8, res8, err8); fprintf('\n');
fprintf(fid, '50,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k9, j9, (k9-1)*50+j9, flag9, ti9, res9, err9); fprintf('\n');
fprintf(fid, '50,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k10(1),k10(2), (k10(1)-1)*50+k10(2), flag10, ti10, res10, err10); fprintf('\n');
fprintf(fid, '60,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k11, j11, (k11-1)*60+j11, flag11, ti11, res11, err11); fprintf('\n');
fprintf(fid, '60,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', k12(1),k12(2), (k12(1)-1)*60+k12(2), flag12, ti12, res12, err12); fprintf('\n');
fprintf(fid, ' 1,  %4i    %4i    %4i   %4i   %8.4f   %10.4e  %10.4e  %11.4e\n', 1,k0, k0, flag0, ti0, res0, err0); fprintf('\n');
%可视化
%t1=1:10:length(resvec1);  
t3=1:10:length(resvec3);
t5=1:10:length(resvec5);  t7=1:10:length(resvec7);
t9=1:10:length(resvec9);  t11=1:10:length(resvec11);
t0=1:10:length(resvec0);
%semilogy(t1,resvec1(t1),'.-k'); 
semilogy(t3,resvec3(t3),'x-k');  hold on
semilogy(t5,resvec5(t5),'<-k');  
semilogy(t7,resvec7(t7),'^-k');  
semilogy(t9,resvec9(t9),'v-k');  
semilogy(t11,resvec11(t11),'>-k'); 
semilogy(t0,resvec0(t0),'*-k');hold off

legend('GMRES(20)','GMRES(30)','GMRES(40)','GMRES(50)','GMRES(60)','  GMRES')
xlabel('迭代次数'); ylabel('相对残差范数');
%axis([1,k0+10,0.5e-11, 1])


 

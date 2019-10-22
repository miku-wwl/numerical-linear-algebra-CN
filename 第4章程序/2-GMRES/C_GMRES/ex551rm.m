clear all; %clc
n=1000;
e=[1:n]'; e1=ones(n-1,1); 
A=diag(e)-diag(e1,-1)-diag(e1,1); A(1,n)=n; A(n,1)=-n;
A=0.5*(A+i*A);
% load SiNa;
% A=Problem.A;  n=size(A,1); A=A+2*eye(n,n);
xs=ones(n,1); b=A*xs;  x=zeros(n,1); tol=1e-12;
max_it=1000;    rm=norm(b); restrt=60;
[x1, k1, j1  time1, res1, resvec1, flag1] = cgmresm (A, b, x, restrt, max_it, tol);       %GMRES,
tic
[x2, flag2, res2, k2, resvec2] = gmres(A, b, restrt, tol, max_it);   %系统自带的GMRES函数
time2=toc;
%显示误差
err1=norm(xs-x1); err2=norm(xs-x2);   
%格式化显示
fprintf('\n');  fid = 1;
fprintf(fid, '   算 法           Out    Int    Iter     flag       Time             Res              Err\n' );     fprintf('\n');
fprintf(fid, 'GMRES0,     %4i    %4i    %4i     %2i     %8.4f      %11.4e  %11.4e\n', k1, j1, (k1-1)*restrt+j1, flag1, time1, res1, err1); fprintf('\n');
fprintf(fid, 'GMRES1,     %4i    %4i    %4i     %2i     %8.4f      %11.4e  %11.4e\n', k2(1),k2(2), (k2(1)-1)*restrt+k2(2), flag2, time2, res2, err2); fprintf('\n');
%可视化
t1=1:length(resvec1);  t2=1:length(resvec2);
semilogy(t1,resvec1,'*-r'); hold on
semilogy(t2,resvec2/rm,'o-'); hold off
legend('GMRES0','GMRES1')


 

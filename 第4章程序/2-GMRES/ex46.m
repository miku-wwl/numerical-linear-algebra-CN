%例4.6-ex46.m
clear all;  
n=1000; e=[1:n]'; e1=ones(n-1,1); 
A=diag(e)+diag(e1,-1)-diag(e1,1); A(1,n)=n; A(n,1)=-n;
% load SiNa; A=Problem.A;  n=size(A,1);  
M=diag(diag(A)); %预处理子
xs=ones(n,1); b=A*xs; x=zeros(n,1); tol=1e-10;
max_it=200;  rm=norm(b);
[x1, k1,  tt1, res1, resvec1, flag1] = mgmres (A, b, x, max_it, tol );        %GMRES,
[x2, k2,  tt2, res2, resvec2, flag2] = pgmres (A, b, x, M, max_it, tol );    %PGMRES,
tic, [x3, flag3, res3, k3, resvec3] = gmres(A, b, [ ], tol, max_it, M); tt3=toc; %系统自带的GMRES函数
%显示误差
err1=norm(xs-x1); err2=norm(xs-x2);  err3=norm(xs-x3); 
%格式化显示
fprintf('\n');  fid = 1;
fprintf(fid, '  算 法          Iter   Iner   flag    Time           Res            Err\n' );     fprintf('\n');
fprintf(fid, ' GMRES,    %4i   %4i    %2i   %8.4f    %11.4e  %11.4e\n', 1, k1, flag1, tt1, res1, err1); fprintf('\n');
fprintf(fid, 'PGMRES,    %4i   %4i    %2i   %8.4f    %11.4e  %11.4e\n',1, k2, flag2, tt2, res2, err2); fprintf('\n');
fprintf(fid, 'XGMRES,    %4i   %4i    %2i   %8.4f    %11.4e  %11.4e\n', k3(1),k3(2), flag2, tt3, res3, err3); fprintf('\n');
%可视化
t1=1:length(resvec1);  t2=1:length(resvec2);  t2=1:length(resvec2);
semilogy(t1,resvec1,'.-k'); hold on
semilogy(t2,resvec2,'*-k'); hold off
%semilogy(t3,resvec3/rm,'*-k'); hold off
legend('GMRES 方法','PGMRES 方法')
xlabel('迭代次数'); ylabel('相对残差范数');
%axis([1,k1+10,0.5e-11, 1])


 

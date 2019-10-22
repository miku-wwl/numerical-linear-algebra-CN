%例4.4-ex44.m
clear all;  
n=1000; e=[1:n]'; e1=ones(n-1,1); 
A=diag(e)+diag(e1,-1)-diag(e1,1); A(1,n)=n; A(n,1)=-n;
% load SiNa; A=Problem.A;  n=size(A,1);  
xs=ones(n,1); b=A*xs; x=zeros(n,1); tol=1e-10;
max_it=200;    rm=norm(b);
[x1, k1,  tt1, res1, resvec1, flag1] = mgmres (A, b, x, max_it, tol );       %GMRES,
tic
[x2, flag2, res2, k2, resvec2] = gmres(A, b, [ ], tol, max_it);   %系统自带的GMRES函数
tt2=toc;
%显示误差
err1=norm(xs-x1); err2=norm(xs-x2);   
%格式化显示
fprintf('\n');  fid = 1;
fprintf(fid, '  算 法          Iter   Iner   flag    Time           Res            Err\n' );     fprintf('\n');
fprintf(fid, 'GMRES1,    %4i   %4i    %2i   %8.4f    %11.4e  %11.4e\n', 1, k1, flag1, tt1, res1, err1); fprintf('\n');
fprintf(fid, 'GMRES2,    %4i   %4i    %2i   %8.4f    %11.4e  %11.4e\n', k2(1),k2(2), flag2, tt2, res2, err2); fprintf('\n');
%可视化
t1=1:length(resvec1);  t2=1:length(resvec2);
semilogy(t1,resvec1,'.-k'); hold on
%semilogy(t2,resvec2/rm,'*-k'); hold off
legend('GMRES 方法')
xlabel('迭代次数'); ylabel('相对残差范数');
axis([1,k1+10,0.5e-11, 1])


 

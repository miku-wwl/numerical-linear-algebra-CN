%例4.3.1
clear all; tol=1.e-10; max_it=500;  
load SiNa; A=Problem.A; n=size(A,1);
xs=ones(n,1); b=A*xs; x0=zeros(n,1); 
[x1,k1,time1,res1,resvec1] = mminres(A,b,x0,max_it,tol);
[x2,k2,time2,res2,resvec2] = minresr(A,b,x0,max_it,tol); %递推形式的MINRES
tic, [x3,flag,res3,k3]= minres(A,b,tol,max_it); time3=toc;  %调用系统函数minres
%格式化显示
fprintf('\n');
fid = 1;
fprintf(fid, '   算  法          Iter         CPU         RES      \n' );     fprintf('\n');
fprintf(fid, 'MMINRES      %3i       %8.4f    %11.4e\n', k1, time1, res1); fprintf('\n');
fprintf(fid, ' MINRES1      %3i       %8.4f    %11.4e\n', k2, time2, res2); fprintf('\n');
fprintf(fid, ' MINRES2      %3i       %8.4f    %11.4e\n', k3, time3, res3); fprintf('\n');
t=1:2:length(resvec1);
semilogy(t,resvec1(t),'.-k');
legend('MINRES')
xlabel('迭代次数'); ylabel('相对残差范数');
er=norm(x1-xs)
%axis([1,65,1.e-10, 1.e4])


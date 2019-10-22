%例4.16
clear all
n=1000; max_it=1000;   tol=1e-10;
A=gallery('lotkin',n); A(1,:)=ones(1,n);
xs=ones(n,1); b=A*xs;   
x0=zeros(n,1); 
[x1,k1,time1,res1,resvec1] = mcgs(A, b, x0, max_it, tol);
tic, [x2,flag,res2,k2,resvec2]=cgs(A, b,  tol, max_it, [ ]); time2=toc; %系统自带
[x,k,time,res,resvec] = mgmres(A, b, x0, max_it, tol);
err1=norm(x1-xs); err2=norm(x2-xs); err=norm(x-xs);     %误差
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '  算 法      Iter    Time      Res          Err\n' );     fprintf('\n');
fprintf(fid, '  mCGS,  %4i     %8.4f   %11.4e    %11.4e\n', k1,  time1, res1, err1); fprintf('\n');
fprintf(fid, '   sCGS,  %4i     %8.4f   %11.4e    %11.4e\n', k2,  time2, res2, err2); fprintf('\n');
fprintf(fid, 'GMRES,  %4i     %8.4f   %11.4e    %11.4e\n', k,  time, res, err); fprintf('\n');
t1=1:length(resvec1);  %t=1:length(resvec);
semilogy(t1,resvec1,'*-k'); %hold on;
%semilogy(t,resvec,'o-k'); hold off
legend('CGS'); xlabel('迭代次数'); ylabel('相对残差范数');
rs1=norm(b-A*x1), rs2=norm(b-A*x2), rs=norm(x-xs),



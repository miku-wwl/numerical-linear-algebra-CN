%例4.15
clear all
n=1000; max_it=1000;   tol=1e-10;
A=gallery('lotkin',n); A(1,:)=ones(1,n);
xs=ones(n,1); b=A*xs; x0=zeros(n,1); 
[x1,k1,time1,res1, resvec1]=bcg(A, b, x0, max_it, tol);   
tic, [x2,flag,res2,k2,resvec2]=bicg(A, b, tol, max_it, [ ]); time2=toc; %系统自带
[x,k,time,res, resvec]=mgmres(A, b, x0, max_it, tol);   
err1=norm(x1-xs); err2=norm(x2-xs); err=norm(x-xs); 
%ILU分解
% setup.type = 'crout';
% setup.milu = 'row';
% setup.droptol = 0.1;
% tic,[M1,M2]=ilu(sparse(A),setup); t3=toc    %ILU分解CPU
% [x3,k3,time3, res3, resvec3] = pbcg(A,b,x0,M1,M2,max_it,tol);
% 
% tic,[x4,flag,res4,k4,resvec4] =  bicg(A, b,  tol, max_it, M1,M2); time4=toc; %系统自带

%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '  算法      Iter   Time        Res        Err\n' );     fprintf('\n');
fprintf(fid, '  BCG,    %4i   %8.4f      %11.4e    %11.4e\n', k1,  time1, res1, err1); fprintf('\n');
fprintf(fid, ' sBiCG,    %4i   %8.4f      %11.4e    %11.4e\n', k2,  time2, res2, err2); fprintf('\n');
fprintf(fid, 'GMRES,   %4i   %8.4f      %11.4e    %11.4e\n', k,  time, res, err); fprintf('\n');
% fprintf(fid, 'sPBiCG,     %4i          %8.4f      %11.4e\n', k4, time4,  res4); fprintf('\n');
t1=1:length(resvec1); 
semilogy(t1,resvec1(t1),'*-k');   hold on
t=1:length(resvec); 
semilogy(t,resvec,'o-k');  hold off
legend('   BCG','GMRES');
xlabel('迭代次数'); ylabel('相对残差范数');
rs1=norm(b-A*x1), rs2=norm(b-A*x2), rs=norm(b-A*x)




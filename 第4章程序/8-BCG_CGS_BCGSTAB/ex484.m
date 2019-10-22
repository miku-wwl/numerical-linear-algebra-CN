%例4.8.4
clear all
n=1000; max_it=1000;   tol=1e-10;
A=gallery('lotkin',n); A(1,:)=ones(1,n);
xs=ones(n,1); b=A*xs;  x0=zeros(n,1); 
setup.type = 'crout';
setup.milu = 'row';
setup.droptol = 0.1;
tic, [M1,M2]=ilu(sparse(A),setup); toc    %ILU分解
[x1,k1,time1,res1,resvec1] = pbcgstab(A, b, x0, M1,M2,max_it, tol);
tic, [x2,flag,res2,k2,resvec2] =bicgstab(A,b,tol,max_it,M1,M2); time2=toc;
err1=norm(x1-xs); err2=norm(x2-xs);   %误差
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '    算   法       Iter       Time              Res                Err\n' );     fprintf('\n');
fprintf(fid, 'PBCGSTAB,  %4i     %8.4f       %11.4e     %11.4e\n', k1,  time1, res1, err1); fprintf('\n');
fprintf(fid, 'sBiCGSTAB,  %4.1f    %8.4f       %11.4e     %11.4e\n', k2, time2,  res2, err2); fprintf('\n');
rs1=norm(b-A*x1), rs2=norm(b-A*x2),   %残差


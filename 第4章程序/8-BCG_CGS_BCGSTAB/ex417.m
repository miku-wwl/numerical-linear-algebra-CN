%例4.17
clear all
n=1000; max_it=1000;   tol=1e-10;
A=gallery('lotkin',n); A(1,:)=ones(1,n);
xs=ones(n,1); b=A*xs;  %M1=diag(diag(A)); M2=tril(A); M3=triu(A); 
x0=zeros(n,1); 
[x1,k1,time1,res1, resvec1]=bcgstab(A, b, x0, max_it, tol);
tic, [x2,flag,res2,k2,resvec2]=bicgstab(A,b,tol,max_it,[ ]); time2=toc; %系统自带
err1=norm(x1-xs); err2=norm(x2-xs);    %误差
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '    算   法     Iter     Time              Res               Err\n' );     fprintf('\n');
fprintf(fid, ' BCGSTAB,  %4i  %8.4f      %11.4e      %11.4e\n', k1,  time1, res1, err1); fprintf('\n');
fprintf(fid, 'sBCGSTAB, %4.1f %8.4f      %11.4e      %11.4e\n', k2,  time2, res2, err2); fprintf('\n');

t1=1:length(resvec1); %t2=1:length(resvec2);
semilogy(t1,resvec1,'*-k');  %hold on
%semilogy(t2,resvec2/norm(b),'+-b'); hold off
legend('BiCGSTAB');
xlabel('迭代次数'); ylabel('相对残差范数');
rs1=norm(b-A*x1),  rs2=norm(b-A*x2),   %残差


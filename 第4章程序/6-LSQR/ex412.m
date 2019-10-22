%例4.12
clear all
n=1000; max_it=1000;   tol=1e-10;
% A=gallery('lotkin',n); A(1,:)=ones(1,n);
e1=ones(n,1);
A=spdiags([-2*e1, 4*e1, -e1], [-1,0,1],n,n);
xs=ones(n,1); b=A*xs;  %M1=diag(diag(A)); M2=tril(A); M3=triu(A); 
x0=zeros(n,1); 
[x1,k1,time1,res1, resvec1] = mlsqr(A, b, x0, max_it, tol);
[x2,k2,time2,res2, resvec2] = mgmres(A, b, x0, max_it, tol);
tic
[x3,f3,res3,k3, resvec3] = lsqr(A, b, tol, max_it, [ ]);  %系统自带
time3=toc;
err1=norm(x1-xs); err2=norm(x2-xs); err3=norm(x3-xs);
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, ' 算   法    Iter        Time            Res                Err\n' );     fprintf('\n');
fprintf(fid, 'mLSQR,  %4i     %8.4f     %11.4e   %11.4e\n', k1,  time1, res1, err1); fprintf('\n');
fprintf(fid, 'GMRES,  %4i     %8.4f     %11.4e    %11.4e\n', k2, time2, res2, err2); fprintf('\n');
fprintf(fid, ' sLSQR,  %4i     %8.4f      %11.4e   %11.4e\n',  k3, time3, res3, err3); fprintf('\n');
t1=1:2:length(resvec1); t2=1:2:length(resvec2); %t3=1:2:length(resvec3);
semilogy(t1,resvec1(t1),'*-k'); hold on
semilogy(t2,resvec2(t2),'o-k');  hold off
%semilogy(t3,resvec3(t3)/norm(b),'p-k');  hold off
legend('LSQR','GMRES');
%legend('mLSQR','GMRES','sLSQR');
xlabel('迭代次数'); ylabel('相对残差范数');
rs1=norm(b-A*x1), rs2=norm(b-A*x2),



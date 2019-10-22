%例4.13
n=1000; max_it=200;   tol=1e-10; e=ones(n,1);
A=spdiags([-e,-2*e,-3*e,12*e,3*e,2*e,e],[-3,-2,-1,0,1,2,3],n,n);
xs=ones(n,1); b=A*xs; x0=zeros(n,1); 
[x1,k1,time1,res1,resvec1] =  gcr(A, b, x0, max_it, tol);
[x2,k2,time2,res2,resvec2] = mgmres(A, b, x0, max_it, tol);
[x3,k3,time3,res3,resvec3] = mlsqr(A, b, x0, max_it, tol);
[x4,k4,time4,res4,resvec4] = mqmres(A, b, x0, max_it, tol);
err1=norm(x1-xs); err2=norm(x2-xs); err3=norm(x3-xs); err4=norm(x4-xs); 
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '  算 法     Iter    Time          Res                Err\n' );     fprintf('\n');
fprintf(fid, '    GCR, %4i  %8.4f   %11.4e    %11.4e\n', k1,time1,res1,err1); fprintf('\n');
fprintf(fid, 'GMRES, %4i  %8.4f   %11.4e    %11.4e\n', k2,time2, res2,err2); fprintf('\n');
fprintf(fid, '  LSQR,  %4i  %8.4f   %11.4e    %11.4e\n', k3,time3, res3,err3); fprintf('\n');
fprintf(fid, 'QMRES, %4i  %8.4f   %11.4e    %11.4e\n', k4,time4, res4,err4); fprintf('\n');

%norm(x)
t1=1:length(resvec1); t2=1:length(resvec2);  t3=1:length(resvec3); t4=1:length(resvec4);
semilogy(t1,resvec1,'+-k');  hold on
semilogy(t2,resvec2,'o-k');   %hold off
semilogy(t3,resvec3,'*-k');   %hold off
semilogy(t4,resvec4,'p-k'); hold off
legend('GCR','GMRES','LSQR','QMRES');
xlabel('迭代次数'); ylabel('相对残差范数');
rs1=norm(b-A*x1), rs2=norm(b-A*x2), rs3=norm(b-A*x3), rs4=norm(b-A*x4),
%axis([1,65,1.e-10, 1.e4])


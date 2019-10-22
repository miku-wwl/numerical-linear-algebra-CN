%例4.14
clear all
n=1000; max_it=200;   tol=1e-10; e=ones(n,1);
A=spdiags([-e,-2*e,-3*e,12*e,3*e,2*e,e],[-3,-2,-1,0,1,2,3],n,n);
xs=ones(n,1); b=A*xs;   
x0=zeros(n,1); restrt=6;
[x1,k1, i1, time1,res1, resvec1] = gcrm(A, b, x0, restrt, max_it, tol);
[x2,k2, i2, time2,res2,resvec2] = gmresm(A, b, x0, restrt, max_it, tol);
err1=norm(x1-xs); err2=norm(x2-xs);
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '  算  法      Out    Int    Iter     Time           Res     Err\n' );     fprintf('\n');
fprintf(fid, ' GCR(6),    %4i   %4i   %4i   %8.4f       %11.4e     %11.4e\n', k1,i1, (k1-1)*restrt+i1, time1, res1, err1); fprintf('\n');
fprintf(fid, 'GMRES(6), %4i   %4i   %4i   %8.4f       %11.4e     %11.4e\n', k2, i2, (k2-1)*restrt+i2, time2, res2, err2); fprintf('\n');
%norm(x)
t1=1:length(resvec1); t2=1:length(resvec2);
semilogy(t1,resvec1,'*-k'); hold on
semilogy(t2,resvec2,'o-k'); hold off
legend('GCR(6)','GMRES(6)');
xlabel('迭代次数'); ylabel('相对残差范数');
rs1=norm(b-A*x1), rs2=norm(b-A*x2),
%axis([1,24,1.e-12, 1.e0])


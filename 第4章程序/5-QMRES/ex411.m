%例4.11
n=1000; max_it=1000;  tol=1e-10;
e=ones(n,1); e1=ones(n-1,1); 
A=diag(4*e)-diag(2*e1,-1)-diag(e1,1); 
%A=diag([1:n]')+diag(e1,-1)-diag(e1,1);  A(1,n)=n; A(n,1)=-n;  %依赖于系数矩阵
xs=e; b=A*xs;  M=tril(A);  M1=diag(diag(A)); 
x0=zeros(n,1); 
[x1,k1,time1,res1,resvec1] = mqmres(A, b, x0, max_it, tol);
[x2,k2,time2,res2,resvec2] = pqmres(A, b, x0, M,max_it, tol);
[x3,k3,time3,res3,resvec3] = mgmres(A, b, x0, max_it, tol);
err1=norm(x1-xs); err2=norm(x2-xs); err3=norm(x3-xs);
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '   算   法      Iter       Time        Res           Err\n' );     fprintf('\n');
fprintf(fid, '  QMRES,  %4i      %8.4f    %11.4e  %11.4e\n', k1,  time1, res1,err1); fprintf('\n');
fprintf(fid, 'PQMRES,  %4i      %8.4f    %11.4e  %11.4e\n', k2, time2, res2,err2); fprintf('\n');
fprintf(fid, '  GMRES,  %4i      %8.4f    %11.4e  %11.4e\n', k3, time3, res3,err3); fprintf('\n');
t1=1:length(resvec1); t2=1:length(resvec2); t3=1:length(resvec3);
semilogy(t1,resvec1,'*-k'); hold on
semilogy(t2,resvec2,'o-k');  
semilogy(t3,resvec3,'p-k'); hold off
legend('  QMRES','PQMRES','  GMRES');
xlabel('迭代次数'); ylabel('相对残差范数');
rs1=norm(b-A*x1), rs2=norm(b-A*x2), rs3=norm(b-A*x3),
%axis([1,65,1.e-10, 1.e4])


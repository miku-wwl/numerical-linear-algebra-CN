%例4.3
n=1000; max_it=500;   tol=1e-10; e=ones(n,1);
A=spdiags([-2*e,-3*e,12*e,3*e,2*e],[-2,-1,0,1,2],n,n);
xs=ones(n,1); b=A*xs;   x0=zeros(n,1); 
[x1,k1,time1,res1,resvec1] = cgnr(A, b, x0, max_it, tol);
[x2,k2,time2,res2,resvec2] = cgne(A, b, x0, max_it, tol);
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '  算   法        Iter           Time               Res\n' );     fprintf('\n');
fprintf(fid, ' CGNR,       %4i         %8.4f      %11.4e\n', k1,  time1,res1); fprintf('\n');
fprintf(fid, ' CGNE,       %4i         %8.4f      %11.4e\n', k2,  time2,  res2); fprintf('\n');
t1=1:length(resvec1); t2=1:length(resvec2);  
semilogy(t1,resvec1,'*-k');  hold on
semilogy(t2,resvec2,'o-k');   hold off
legend('CGNR','CGNE');
xlabel('迭代次数'); ylabel('相对残差范数');
er1=norm(x1-xs),   er2=norm(x2-xs),
rs1=norm(b-A*x1), rs2=norm(b-A*x2), 
%axis([1,65,1.e-10, 1.e4])


%例4.12b
clear all
n=1000; max_it=1000;   tol=1e-10;
% A=gallery('lotkin',n); A(1,:)=ones(1,n);
e1=ones(n,1);
A=spdiags([-2*e1, 4*e1, -e1], [-1,0,1],n,n);
M=diag(diag(A))+diag(diag(A,-1),-1);
xs=ones(n,1); b=A*xs;  %M1=diag(diag(A)); M2=tril(A); M3=triu(A); 
x0=zeros(n,1);    rm=norm(b);
tic
[x3,f3,res3,k3, resvec3] = lsqr(A, b, tol, max_it, M);  %系统自带 A*inv(M)*y=b; y=M*x
time3=toc;
%[x1,k1,time1,res1, resvec1] = plsqr(A, b, x0, M,max_it, tol);
[x1,k1,time1,res1, resvec1] = plsqr(A, b, x0, M,max_it, tol);
[x2,k2,time2,res2, resvec2] = mlsqr(A, b, x0, max_it, tol);
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '   算   法        Iter           Time                   Err\n' );     fprintf('\n');
fprintf(fid, 'sLSQR,     %4i         %8.4f      %11.4e\n', k3, time3, res3); fprintf('\n');
fprintf(fid, 'PLSQR,     %4i        %8.4f       %11.4e\n', k1,  time1, res1); fprintf('\n');
fprintf(fid, '  LSQR,     %4i         %8.4f      %11.4e\n', k2, time2, res2); fprintf('\n');
%norm(x)
t1=1:length(resvec1); t2=1:length(resvec2);
t3=1:length(resvec3);
semilogy(t1,resvec1,'*-r'); hold on
semilogy(t2,resvec2,'o-b'); 
semilogy(t3,resvec3/rm,'d-g'); hold off
legend('sLSQR','LSQR','PLSQR');
rs1=norm(b-A*x1)
er1=norm(x1-xs)
rs2=norm(b-A*x2)
er2=norm(x2-xs)
%axis([1,65,1.e-10, 1.e4])


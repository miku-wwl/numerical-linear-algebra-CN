%例4.10
% n=1000; tol=1e-10; max_it=100;
% A=zeros(n,n);
% for i=1:n
%     for j=1:n
%         A(i,j)=0.5/(n-i-j+1.5);
%     end
% end
clear all
tol=1.e-10; max_it=1000;  
load SiNa;
A=Problem.A;  n=size(A,1);
xs=ones(n,1); b=A*xs; 
x0=zeros(n,1); 
[x1,k1,time1,res1,resvec1] = msymmlq(A,b,x0,max_it,tol);
[x2,k2,time2,res2,resvec2] = mminres(A, b, x0,max_it, tol);
[x3,k3,time3,res3,resvec3] = mgmres(A, b, x0,max_it, tol);
tic; [x4,flag,res4,k4,resvec4] = symmlq(A,b,tol,max_it,[],[],x0); time4=toc; %[x,flag,relres,iter,resvec] = symmlq(A,b,tol,maxit,M1,M2,x0)
err1=norm(x1-xs);  err2=norm(x2-xs);  err3=norm(x3-xs);  err4=norm(x4-xs);
%格式化显示
fprintf('\n'); fid = 1;
fprintf(fid, '   算   法           Iter           Time         Res          Err\n' );     fprintf('\n');
fprintf(fid, 'SYMMLQ,    %4i         %8.4f      %11.4e   %11.4e\n', k1,  time1, res1, err1); fprintf('\n');
fprintf(fid, ' MINRES,     %4i         %8.4f      %11.4e   %11.4e\n', k2, time2,  res2, err2); fprintf('\n');
fprintf(fid, ' GMRES,      %4i         %8.4f      %11.4e   %11.4e\n', k3, time3,  res3, err3); fprintf('\n');
fprintf(fid, ' XSYMLQ,    %4i         %8.4f      %11.4e   %11.4e\n', k4, time4,  res4, err4); fprintf('\n');
%norm(x)
t1=1:10:length(resvec1); t2=1:10:length(resvec2); t3=1:10:length(resvec3); 
semilogy(t1,resvec1(t1),'o-r');  hold on
semilogy(t2,resvec1(t2),'+-b');  
semilogy(t3,resvec3(t3),'*-g');   
%beta=norm(b); resvec4=resvec4/beta; t4=1:10:length(resvec4);
%semilogy(t4,resvec4(t4),'v-k'); hold off
legend('SYMMLQ','MINRES','GMRES');
rs1=norm(b-A*x1)
%axis([1,65,1.e-10, 1.e4])


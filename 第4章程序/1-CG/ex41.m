%例4.1
clear all
n=1000; tol=1e-10;  max_it=1000;
e=[1:n]'; e1=ones(n-1,1);
A=diag(e)-diag(e1,1)-diag(e1,-1);
xs=ones(n,1);  b=A*xs;  x0=zeros(n,1);
[x, k ,time, res, resvec] = mcg(A,b,x0,max_it,tol);
%格式化显示
fprintf('\n');
fid = 1;
fprintf(fid, '算 法     Iter         CPU         RES      \n' );     fprintf('\n');
fprintf(fid, 'CG      %3i       %8.4f    %11.4e\n', k, time, res); fprintf('\n');
%可视化
t=1:length(resvec);
semilogy(t,resvec,'*-r'); 
legend('CG 方法'); axis([1,k,1.e-12, 1.e2])
xlabel('迭代次数'); ylabel('相对残差范数');
err=norm(x-xs)



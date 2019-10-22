%例4.3.2
clear all
tol=1.e-10; max_it=500;  
load SiNa; A=Problem.A;  n=size(A,1); %n=5743
mu=-1.0;
tic
L=ichol(sparse(A-mu*eye(n)));  
M1=L; M2=L';   M=M1*M2;
toc
%M1=eye(n);
%M3=diag(diag(A));  %Jacobi预条件子对此例无效
%D=diag(diag(A)); L=-tril(A,-1); U=triu(A,1); w=1.05;
%M4=1/(w*(2-w))*(D-w*L)*(D\(D-w*U)); %SSOR预条件子对本例无效
%[V,D]=eig(A); M3=V*V';  %特征分解对本例无效 
xs=ones(n,1); b=A*xs; 
x0=zeros(n,1); 
[x1,k1,time1,r1,rvec1] = mminres(A,b,x0,max_it,tol);    
%[x2,k2,time2,r2,rvec2] = pminres(A,b,x0,M,max_it,tol);  
[x3,k3,time3,r3,rvec3] = pminres1(A,b,x0,M1,M2,max_it,tol);
tic, [x4,flag,r4,k4]= minres(A,b,tol,max_it,M1,M2); time4=toc;  %调用系统函数minres
%可视化
t1=1:5:length(rvec1);
semilogy(t1,rvec1(t1),'-^b'); hold on
t3=1:2:length(rvec3);
semilogy(t3,rvec3(t3),'-*r'); hold off
legend('MINRES','PMINRES')
xlabel('迭代次数'); ylabel('相对残差范数');
axis([1, k1, 1.e-11, 1])
%格式化显示
fprintf('\n');
fid = 1;
fprintf(fid, '   算   法           Iter          CPU          RES      \n' );     fprintf('\n');
fprintf(fid, 'MMINRES      %3i       %8.4f    %11.4e\n', k1, time1, r1); fprintf('\n');
%fprintf(fid, ' PMINRES      %3i       %8.4f    %11.4e\n', k2, time2, r2); fprintf('\n');
fprintf(fid, 'PMINRES1      %3i       %8.4f    %11.4e\n', k3, time3, r3); fprintf('\n');
fprintf(fid, 'PMINRES2      %3i       %8.4f    %11.4e\n', k4, time4, r4); fprintf('\n');
%误差
er1=norm(x1-xs)
er2=norm(x3-xs)

% A=(L\A)/L'; b=L\b; [y,k,time,r,rvec]=mminres(A,b,x0,max_it,tol); x=L'\y; k,time,  

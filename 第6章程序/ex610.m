%例6.10-ex610.m
clear all; n=1024; 
e=ones(n,1); e1=ones(1,n-1);
b=4*e; a=[0, -2*e1]'; c=[-e1,0]';
f=e; f(1)=3; f(n)=2;  l1=1; u1=1;
tic, [x1]=mchase(a,b,c,f);  t1=toc;
tic, [x2]=mchase_var(a,b,c,f,l1,u1);  t2=toc;
A=4*diag(e)-2*diag(e1,-1)-diag(e1,1);
r1=norm(A*x1-f);  r2=norm(A*x2-f); 
fprintf('\n'); fid = 1;
fprintf(fid, '   算法       误差     CPU时间   ' );     fprintf('\n');
fprintf(fid, '追 赶  法     %8.4e     %8.4f\n', r1, t1);   
fprintf(fid, '变参数法     %8.4e     %8.4f\n', r2, t2);    fprintf('\n');

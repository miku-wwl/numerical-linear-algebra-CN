%例3.7(1)-ex37_1.m
clear all
m=2000; t=3; n=m*t; tol=1.e-12; a=1;
B=cell(m,1); A=cell(m,1); C=cell(m,1);
for i=1:m
    B{i}=[4,-1,0; -1,4,-1; 0,-1,4];
    A{i}=diag([-1,-1,-1]); C{i}=A{i};
    f{i}=[1, 2, 3]'+(i-1)*t;  %右端向量
    x0{i}=zeros(t,1); %初始向量
end
[x,iter1,res1,time1]=PEa1(A,B,C,f,0.5,x0,tol);
[x,iter2,res2,time2]=PEa1(A,B,C,f,1.0,x0,tol);
[x,iter3,res3,time3]=PEa1(A,B,C,f,1.4,x0,tol);
[x,iter4,res4,time4]=PEa1(A,B,C,f,1.5,x0,tol);
[x,iter5,res5,time5]=PEa1(A,B,C,f,1.6,x0,tol);
[x,iter6,res6,time6]=PEa1(A,B,C,f,2.0,x0,tol);
[x,iter7,res7,time7]=PEa1(A,B,C,f,0.0,x0,tol);
%格式化显示
fid = 1; 
fprintf(fid, ' PE方法         Iter      CPU           RES\n' );     fprintf('\n');
fprintf(fid, ' a=0.5     %4i     %4.4f     %11.4e\n', iter1,time1,res1); fprintf('\n');
fprintf(fid, ' a=1.0     %4i     %4.4f     %11.4e\n', iter2,time2,res2); fprintf('\n');
fprintf(fid, ' a=1.4     %4i     %4.4f     %11.4e\n', iter3,time3,res3); fprintf('\n');
fprintf(fid, ' a=1.5     %4i     %4.4f     %11.4e\n', iter4,time4,res4); fprintf('\n');
fprintf(fid, ' a=1.6     %4i     %4.4f     %11.4e\n', iter5,time5,res5); fprintf('\n');
fprintf(fid, ' a=2.0     %4i     %4.4f     %11.4e\n', iter6,time6,res6); fprintf('\n');
fprintf(fid, ' a=0.0     %4i     %4.4f     %11.4e\n', iter7,time7,res7); fprintf('\n');




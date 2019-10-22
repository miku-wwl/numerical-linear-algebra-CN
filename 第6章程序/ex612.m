%例6.12-ex612.m
clear all; m=1000; 
A=cell(m,1); B=cell(m,1); C=cell(m,1); f=cell(m,1);
Bi=[4,-1,0;-1,4,-1;0,-1,4];
Ai=[13,0,0;0,11,0;1,0,12]; Ci=Ai'; fi=[1,0,1]';
for i=1:m
    A{i}=Ai; B{i}=Bi; C{i}=Ci; f{i}=fi;
end
tic, [x1]=mchase_block(A,B,C,f,m); t1=toc;
tic, [x2]=mdouble_par(A,B,C,f,m);  t2=toc;
res1(1)=norm(B{1}*x1{1}+C{1}*x1{2}-f{1});
res2(1)=norm(B{1}*x2{1}+C{1}*x2{2}-f{1});
for i=2:m-1
    res1(i)=norm(A{i}*x1{i-1}+B{i}*x1{i}+C{i}*x1{i+1}-f{i});
    res2(i)=norm(A{i}*x2{i-1}+B{i}*x2{i}+C{i}*x2{i+1}-f{i});
end
res1(m)=norm(A{m}*x1{m-1}+B{m}*x1{m}-f{m});
res2(m)=norm(A{m}*x2{m-1}+B{m}*x2{m}-f{m});
r1=max(res1); r2=max(res2);
fprintf('\n'); fid = 1;
fprintf(fid, '   算法           误差            CPU时间   ' );     fprintf('\n');
fprintf(fid, '块追赶法     %8.4e   %8.4f\n', r1, t1);   
fprintf(fid, '双参数法     %8.4e   %8.4f\n', r2, t2);    fprintf('\n');

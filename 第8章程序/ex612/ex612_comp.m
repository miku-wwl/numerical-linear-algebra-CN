function ex612_comp
% 将Jacobi,GS,BJ,BGS,QGS五个算法进行对比
clear
n = 60; n=120; n=180; n=240; n=300; n=360;
ni=3; m=n/ni; I=eye(n);
A = cell(m,m);  B = cell(m,m);  Ii = -eye(ni);
AI = [4 -1 -1; -1  4 -1; -1 -1  4];  
BI = [4 -1 0; -1  4 -1; 0 -1  4]; 
tic
for i=1:m
    for j=1:m
        A{i,j}=zeros(ni); B{i,j}=zeros(ni);
    end
end
for i=2:m-1
     A{i,i-1}=Ii; A{i,i}=AI; A{i,i+1}=Ii;
     B{i,i-1}=Ii; B{i,i}=BI; B{i,i+1}=Ii;
end
A{1,1}=AI;  A{1,2}=Ii; A{m,m-1}=Ii; A{m,m}=AI;
B{1,1}=BI;  B{1,2}=Ii; B{m,m-1}=Ii; B{m,m}=BI;
A=cell2mat(A);
B=cell2mat(B);

%n=120;
[Itr1, Err1, t1]=ex612J(A,B);
[Itr2, Err2, t2]=ex612GS(A,B);
[Itr3, Err3, t3]=ex612QGS(A,B) ;
[Itr4, Err4, t4]=ex612BJ(A,B);
[Itr5, Err5, t5]=ex612BGS(A,B);
% 输出结果
fid = 1;
logBody = ' %2i   %8.4e  %8.4f    \n';
fprintf('  Itr            Err           time    '); fprintf('\n');
fprintf(logBody, [Itr1, Err1, t1]); fprintf('\n');
fprintf(logBody, [Itr2, Err2, t2]); fprintf('\n');
fprintf(logBody,[Itr3, Err3, t3]); fprintf('\n');
fprintf(logBody,[Itr4, Err4, t4]); fprintf('\n');
fprintf(logBody,[Itr5, Err5, t5]); fprintf('\n');



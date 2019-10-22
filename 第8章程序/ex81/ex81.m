%Àý8.1
clc;  
ni=3; m=100;   
A = cell(m,m);  B = cell(m,m);  Ii = -eye(ni);
AI = [4 -1 -1;  -1  4 -1; -1 -1  4];   
BI =  [4 -1 0;  -1  4 -1; 0 -1  4]; 
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
A=cell2mat(A); B=cell2mat(B); F=eye(size(A,1));
tic
[X1]=BSM(A,A',F);
toc
tic
[X2]=BSM(A,B,F);
toc
Err1=norm(A*X1+X1*A'-F,'fro')
Err2=norm(A*X2+X2*B-F,'fro')

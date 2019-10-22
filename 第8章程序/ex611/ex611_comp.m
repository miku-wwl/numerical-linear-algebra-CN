function ex611_comp 
%clc;  
alpha=[0.2:0.3:0.8, 0.9:0.1:1.5,  2  4  6 12];
for i=1:14
    [k(i), t(i)]=ex611(alpha(i));
end
%格式化显示
fid = 1;
fprintf(fid, '算法     Iter     CPU \n' );     fprintf('\n');
fprintf(fid, 'alpha(1), %4i    %4.4f \n', k(1),t(1)); fprintf('\n');
fprintf(fid, 'alpha(2), %4i    %4.4f \n', k(2),t(2)); fprintf('\n');
fprintf(fid, 'alpha(3), %4i    %4.4f \n', k(3),t(3)); fprintf('\n');
fprintf(fid, 'alpha(4), %4i    %4.4f \n', k(4),t(4)); fprintf('\n');
fprintf(fid, 'alpha(5), %4i    %4.4f \n', k(5),t(5)); fprintf('\n');
fprintf(fid, 'alpha(6), %4i    %4.4f \n', k(6),t(6)); fprintf('\n');
fprintf(fid, 'alpha(7), %4i    %4.4f \n', k(7),t(7)); fprintf('\n');
fprintf(fid, 'alpha(8), %4i    %4.4f \n', k(8),t(8)); fprintf('\n');
fprintf(fid, 'alpha(9), %4i    %4.4f \n', k(9),t(9)); fprintf('\n');
fprintf(fid, 'alpha(10), %4i    %4.4f \n', k(10),t(10)); fprintf('\n');
fprintf(fid, 'alpha(11), %4i    %4.4f \n', k(11),t(11)); fprintf('\n');
fprintf(fid, 'alpha(12), %4i    %4.4f \n', k(12),t(12)); fprintf('\n');
fprintf(fid, 'alpha(13), %4i    %4.4f \n', k(13),t(13)); fprintf('\n');
fprintf(fid, 'alpha(14), %4i    %4.4f \n', k(14),t(14)); fprintf('\n');

%计算最优参数
ni=3; m=200;  I=eye(ni*m); F=ones(ni*m);
A = cell(m,m);  B = cell(m,m);  Ii = -eye(ni);
C = cell(m,m);  D = cell(m,m); 
AI = [4 -1 0;  -1  4 -1; 0 -1  4];   
BI =  [4 1 1;  1  4 1; 1 1  4]; 
CI=[4.5 -1 -1; -1 4.5 -1; -1 -1 4.5];
DI= [4 -1 -1;  -1  4 -1; -1 -1  4]; 
for i=1:m
    for j=1:m
        A{i,j}=zeros(ni); B{i,j}=zeros(ni);
        C{i,j}=zeros(ni); D{i,j}=zeros(ni);
    end
end
for i=2:m-1
     A{i,i-1}=Ii; A{i,i}=AI; A{i,i+1}=Ii;
     B{i,i-1}=Ii; B{i,i}=BI; B{i,i+1}=Ii;
     C{i,i-1}=Ii; C{i,i}=CI; C{i,i+1}=Ii;
     D{i,i-1}=Ii; D{i,i}=DI; D{i,i+1}=Ii;
end
A{1,1}=AI;  A{1,2}=Ii; A{m,m-1}=Ii; A{m,m}=AI;
B{1,1}=BI;  B{1,2}=Ii; B{m,m-1}=Ii; B{m,m}=BI;
C{1,1}=CI;  C{1,2}=Ii; C{m,m-1}=Ii; C{m,m}=CI;
D{1,1}=DI;  D{1,2}=Ii; D{m,m-1}=Ii; D{m,m}=DI;
A=cell2mat(A); B=cell2mat(B);
C=cell2mat(C); D=cell2mat(D);
eta=eig(inv(A)*C);
mu=eig(B*inv(D));
ap=sqrt(max(eta)*min(eta))
aq=sqrt(max(mu)*min(mu))
b=eig(B);
d=eig(D);
a=eig(A);
c=eig(C);
if ap<=aq,
    apw=sqrt(max(c)*min(c)/(max(a)*min(a)))
else
    aqw=sqrt(max(b)*min(b)/(max(d)*min(d)))
end

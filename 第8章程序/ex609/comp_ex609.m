%clc
%format short e
alpha=[0.2:0.1:1.0, 2:2:10,15,20];
for i=1:16
    [k1(i), t1(i),k2(i),t2(i)]=ex609(alpha(i));
end
%格式化显示
fid = 1;
fprintf(fid, '算法     Iter     CPU \n' );     fprintf('\n');
fprintf(fid, 'alpha(1), %4i    %4.4f , %4i    %4.4f \n',k1(1),t1(1), k2(1),t2(1)); fprintf('\n');
fprintf(fid, 'alpha(2), %4i    %4.4f , %4i    %4.4f\n', k1(2),t1(2), k2(2),t2(2)); fprintf('\n');
fprintf(fid, 'alpha(3), %4i    %4.4f , %4i    %4.4f\n', k1(3),t1(3), k2(3),t2(3)); fprintf('\n');
fprintf(fid, 'alpha(4), %4i    %4.4f , %4i    %4.4f\n', k1(4),t1(4), k2(4),t2(4)); fprintf('\n');
fprintf(fid, 'alpha(5), %4i    %4.4f, %4i    %4.4f \n', k1(5),t1(5), k2(5),t2(5)); fprintf('\n');
fprintf(fid, 'alpha(6), %4i    %4.4f , %4i    %4.4f\n', k1(6),t1(6), k2(6),t2(6)); fprintf('\n');
fprintf(fid, 'alpha(7), %4i    %4.4f, %4i    %4.4f \n', k1(7),t1(7), k2(7),t2(7)); fprintf('\n');
fprintf(fid, 'alpha(8), %4i    %4.4f, %4i    %4.4f \n', k1(8),t1(8), k2(8),t2(8)); fprintf('\n');
fprintf(fid, 'alpha(9), %4i    %4.4f , %4i    %4.4f\n', k1(9),t1(9), k2(9),t2(9)); fprintf('\n');
fprintf(fid, 'alpha(10), %4i    %4.4f, %4i    %4.4f \n', k1(10),t1(10),  k2(10),t2(10)); fprintf('\n');
fprintf(fid, 'alpha(11), %4i    %4.4f , %4i    %4.4f\n',k1(11),t1(11),  k2(11),t2(11)); fprintf('\n');
fprintf(fid, 'alpha(12), %4i    %4.4f , %4i    %4.4f\n', k1(12),t1(12), k2(12),t2(12)); fprintf('\n');
fprintf(fid, 'alpha(13), %4i    %4.4f, %4i    %4.4f \n', k1(13),t1(13), k2(11),t2(11)); fprintf('\n');
fprintf(fid, 'alpha(14), %4i    %4.4f, %4i    %4.4f \n', k1(14),t1(14), k2(12),t2(12)); fprintf('\n');
fprintf(fid, 'alpha(15), %4i    %4.4f, %4i    %4.4f \n', k1(15),t1(15), k2(11),t2(11)); fprintf('\n');
fprintf(fid, 'alpha(16), %4i    %4.4f, %4i    %4.4f \n', k1(16),t1(16), k2(12),t2(12)); fprintf('\n');
%计算近似最优参数
ni=3; m=100;  I=eye(ni*m); F=I;
A = cell(m,m);  B = cell(m,m);  Ii = -eye(ni);
AI = [12 1 1;  1  17 0; 1 0  15];  
BI = [9 1 1; 1 8 1; 1 1 7]; 
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
lam2=min(eig(A)); lam1=max(eig(A));
mu2=min(eig(B)); mu1=max(eig(B));
aopt=sqrt((mu1*mu2)/(lam1*lam2))


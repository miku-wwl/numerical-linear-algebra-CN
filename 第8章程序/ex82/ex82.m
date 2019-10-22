clear all;  %调用时要注意alpha是负数
ni=3; m=100;  I=eye(ni*m); F=I;
A = cell(m,m);  B = cell(m,m);  Ii = -eye(ni);
Ai = [4 -1 -1; -1  4 -1; -1 -1  4];  
Bi = [4 -1 0; -1  4 -1; 0 -1  4]; 
for i=1:m
    for j=1:m
        A{i,j}=zeros(ni); B{i,j}=zeros(ni);
    end
end
for i=2:m-1
     A{i,i-1}=Ii; A{i,i}=Ai; A{i,i+1}=Ii;
     B{i,i-1}=Ii; B{i,i}=Bi; B{i,i+1}=Ii;
end
A{1,1}=Ai;  A{1,2}=Ii; A{m,m-1}=Ii; A{m,m}=Ai;
B{1,1}=Bi;  B{1,2}=Ii; B{m,m-1}=Ii; B{m,m}=Bi;
A=cell2mat(A); B=cell2mat(B);
alpha=-[0.1:0.1:1.0, 5,8];
for i=1:12
    [k(i), t(i)]=par_iter(A,B,F,alpha(i));
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
%计算最忧参数      
Ae=eig(A); Be=eig(B);
lam2=min(Ae); lam1=max(Ae);
mu2=min(Be); mu1=max(Be);
alam=-1./sqrt(lam1*lam2)
almu=-1./sqrt(mu1*mu2) 

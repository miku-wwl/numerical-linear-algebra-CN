%例7.8
clear all
A=[3 2 3 4 5 6 7;11 1 2 3 4 5 6; 2 8 9 1 2 3 4;-4 2 9 11 13 15 8; 
    -1 -2 -3 -1 -1 -1 -1; 3 2 3 4 13 15 8; -2 -2 -3 -4 -5 -3 -3];
tic,  [iter1,D1]=qr_eig(A);  t1=toc;
tic,  [iter2,D2]=ddiqr_eig(A);  t2=toc;
fprintf('\n'); fid = 1;
fprintf(fid, '   算    法         迭代次数         CPU时间   ' );     fprintf('\n');
fprintf(fid, ' 基本QR方法     %4i             %8.4f\n', iter1, t1);   
fprintf(fid, '双位移QR方法   %4i               %8.4f\n', iter2,t2);   fprintf('\n');
disp( '  基本QR方法特征值,   双位移QR方法特征值')
disp([D1,       D2])
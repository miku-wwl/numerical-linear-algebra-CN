clc
%format short e
alpha = -[0.1:0.1:1.0, 3, 5, 8];
for i=1:13
    [k1(i), k2(i), t1(i), t2(i)]=ex617(alpha(i));
end
%格式化显示
fid = 1;
fprintf(fid, '算法       Iter1   Iter2    CPU1   CPU2\n' );     fprintf('\n');
fprintf(fid, 'alpha(1), %4i    %4i   %4.4f   %4.4f\n', k1(1),k2(1), t1(1), t2(1)); fprintf('\n');
fprintf(fid, 'alpha(2), %4i    %4i   %4.4f   %4.4f\n', k1(2),k2(2), t1(2), t2(2)); fprintf('\n');
fprintf(fid, 'alpha(3), %4i    %4i  %4.4f    %4.4f\n', k1(3),k2(3), t1(3), t2(3)); fprintf('\n');
fprintf(fid, 'alpha(4), %4i    %4i   %4.4f   %4.4f\n', k1(4),k2(4), t1(4), t2(4)); fprintf('\n');
fprintf(fid, 'alpha(5), %4i    %4i   %4.4f   %4.4f\n', k1(5),k2(5), t1(5), t2(5)); fprintf('\n');
fprintf(fid, 'alpha(6), %4i    %4i   %4.4f   %4.4f\n', k1(6),k2(6), t1(6), t2(6)); fprintf('\n');
fprintf(fid, 'alpha(7), %4i    %4i   %4.4f   %4.4f\n', k1(7),k2(7), t1(7), t2(7)); fprintf('\n');
fprintf(fid, 'alpha(8), %4i    %4i   %4.4f   %4.4f\n', k1(8),k2(8), t1(8), t2(8)); fprintf('\n');
fprintf(fid, 'alpha(9), %4i    %4i   %4.4f   %4.4f\n', k1(9),k2(9), t1(9), t2(9)); fprintf('\n');
fprintf(fid, 'alpha(10), %4i   %4i   %4.4f  %4.4f\n', k1(10),k2(10), t1(10), t2(10)); fprintf('\n');
fprintf(fid, 'alpha(11), %4i   %4i   %4.4f  %4.4f\n', k1(11),k2(11), t1(11), t2(11)); fprintf('\n');
fprintf(fid, 'alpha(12), %4i   %4i   %4.4f   %4.4f\n', k1(12),k2(12), t1(12), t2(12)); fprintf('\n');
fprintf(fid, 'alpha(13), %4i    %4i  %4.4f   %4.4f\n',k1(13),k2(13), t1(13), t2(13)); fprintf('\n');
 

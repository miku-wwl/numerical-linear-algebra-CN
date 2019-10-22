%clc
%format short e
alpha=-[0.1:0.1:1.0, 5,10];
for i=1:12
    [k(i), t(i)]=ex606(alpha(i));
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


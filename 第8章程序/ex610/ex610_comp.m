function ex610_comp 
%clc;  
A=[2 2 1; 2 5 4;1 4 5]; B=[4 1 0; 1 4 1;0 1 3];
C=[2 1 0; 1 2 2; 0 2 3]; D=[2 0 1; 0 2 2; 1 2 3];
F=eye(3);
eta=eig(inv(A)*C)
mu=eig(B*inv(D))
ap=sqrt(max(eta)*min(eta))
aq=sqrt(max(mu)*min(mu))
rhoP=(sqrt(max(eta))-sqrt(min(eta)))/(sqrt(max(eta))+sqrt(min(eta)))
rhoQ=(sqrt(max(mu))-sqrt(min(mu)))/(sqrt(max(mu))+sqrt(min(mu)))
rhoPQp=(ap-min(eta))*(max(mu)-ap)/((ap+min(eta))*(max(mu)+ap))
rhoPQq=(aq-min(eta))*(max(mu)-aq)/((aq+min(eta))*(max(mu)+aq))
%b=eig(B)
%d=eig(D)
a=eig(A)
c=eig(C)
apw=sqrt(max(c)*min(c)/(max(a)*min(a)))
rhoPQpw=(apw-min(eta))*(max(mu)-apw)/((apw+min(eta))*(max(mu)+apw))
%
alpha=[ap, aq, apw, 0.2 0.4 0.6 0.8 1.0 1.5 3.0 6.0  12.0];
for i=1:12
    [k(i), t(i)]=ex610(alpha(i));
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
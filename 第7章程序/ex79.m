%ex7.9
clear all
A=[3 2 3 4 5 6 7; 11 1 2 3 4 5 6; 2 8 9 1 2 3 4;
     -4 2 9 11 13 15 8; -1 -2 -3 -1 -1 -1 -1; 
      3 2 3 4 13 15 8; -2 -2 -3 -4 -5 -3 -3];
[Lam, V, iter, ki]=ddiqr_eigvec(A);
D=eig(A); %用系统函数求A的全部特征值
disp('  双位移QR方法结果   eig函数计算结果')
disp([Lam,  D])  %显示结果
iter  %显示迭代次数
ki   %显示反幂法求每个特征向量的迭代次数
err=norm(A*V-V*diag(Lam),inf)  %验证计算是否正确

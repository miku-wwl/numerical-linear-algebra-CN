%ex438
clear all
n=15; e=[1:n]'; e1=ones(n-1,1); e2=ones(n-2,1); 
e3=ones(n-3,1); e4=ones(n-4,1);
A=diag(e)+diag(e1,1)-diag(e1,-1)-2*diag(e2,2)+2*diag(e2,-2);
A=A+3*diag(e3,3)-3*diag(e3,-3)-2*diag(e4,4)+2*diag(e4,-4);
[Lam, V, iter, ki]=ddiqrm_vec(A);
D=eig(A); %用系统函数求A的全部特征值
disp([Lam,  D])  %显示结果
iter   %显示迭代次数
ki   %显示反幂法求每个特征向量的迭代次数
err=norm(A*V-V*diag(Lam),inf)  %验证计算是否正确
